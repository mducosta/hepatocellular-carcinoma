# ============================================================
# Analise de Expressao Diferencial - LIHC vs. Normal (GTEX)
# Via de lipideos e aterosclerose (KEGG hsa05417)
# ============================================================
#
# OBS.: Esta e a versao ORIGINAL (legado) do estudo.
# Para a versao PORTATIL e REPRODUTIVEL, use o script na raiz:
#     pipeline_hepato.R
#
# ============================================================

# ------------------------------------------------------------------
# 1) Pacotes
# ------------------------------------------------------------------
setwd("C:/Users/oorie/OneDrive/Documentos/TRABALHOS/EXPRESSÃO DIFERENCIAL/LIVER")

suppressPackageStartupMessages({
  library(KEGGREST)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(rio)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(ggraph)
  library(plotly)
  library(scales)
  library(tidyr)
  library(ggalluvial)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Criar pastas de saida
dir.create("results",          showWarnings = FALSE)
dir.create("kegg_enrichment",  showWarnings = FALSE)

# ------------------------------------------------------------------
# 2) Genes da via hsa05417 (KEGG)
# ------------------------------------------------------------------
# shell.exec("https://www.kegg.jp/pathway/hsa05417")
pathway_id <- "hsa05417"

pw <- keggGet(pathway_id)
stopifnot(length(pw) == 1L)

# Usar nome diferente de 'g' para nao conflitar com o igraph adiante
pathway_genes_raw <- pw[[1]]$GENE
if (is.null(pathway_genes_raw) || length(pathway_genes_raw) == 0L)
  stop("No GENE field returned by KEGG for ", pathway_id)

gene_symbols_kegg <- pathway_genes_raw[seq(2, length(pathway_genes_raw), by = 2)] |>
  str_remove("\\s*\\[.*$") |>
  str_trim() |>
  str_extract("^[^;]+") |>
  str_trim() |>
  unique() |>
  sort()

genes_via_tbl <- tibble(
  pathway_id  = pathway_id,
  gene_symbol = gene_symbols_kegg
)

print(genes_via_tbl, n = Inf)
rio::export(genes_via_tbl, "Hsa_genes_lipideos_aterosclerose.csv")

# Vetor de simbolos - usado nas analises seguintes
genes_via <- genes_via_tbl$gene_symbol

# ------------------------------------------------------------------
# 3) Importar tabela do Xena (log2(norm_count + 1))
#    Estrutura esperada:
#      cols 1:5 = metadados (sample, primary_site, sample_type, study,
#                            TCGA_GTEX_main_category)
#      cols 6:n = simbolos genicos
# ------------------------------------------------------------------
liver <- rio::import("liver_lip_aterosclero.tsv")
stopifnot(is.data.frame(liver))

if ("samples" %in% colnames(liver)) liver$samples <- NULL
names(liver) <- sub("^_", "", names(liver))

View(liver)

# ------------------------------------------------------------------
# 4) Verificacoes de consistencia
# ------------------------------------------------------------------
required_cols <- c("sample", "primary_site", "sample_type", "study",
                   "TCGA_GTEX_main_category")
missing_cols  <- setdiff(required_cols, colnames(liver))
if (length(missing_cols) > 0L)
  stop("Missing required columns in `liver`: ", paste(missing_cols, collapse = ", "))

stopifnot(
  all(c("GTEX Liver", "TCGA Liver Hepatocellular Carcinoma") %in%
        unique(liver$TCGA_GTEX_main_category))
)

# ------------------------------------------------------------------
# 5) Separar metadados e expressao; definir condicoes
# ------------------------------------------------------------------
meta <- liver[, required_cols, drop = FALSE]

expr <- liver[, setdiff(colnames(liver), required_cols), drop = FALSE] |>
  as.matrix()
storage.mode(expr) <- "numeric"

meta$condition <- dplyr::case_when(
  meta$TCGA_GTEX_main_category == "GTEX Liver"                          ~ "Normal",
  meta$TCGA_GTEX_main_category == "TCGA Liver Hepatocellular Carcinoma" ~ "LIHC",
  TRUE ~ NA_character_
)
meta$condition <- factor(meta$condition, levels = c("Normal", "LIHC"))

keep_samples <- !is.na(meta$condition)
meta <- meta[keep_samples, , drop = FALSE]
expr <- expr[keep_samples, , drop = FALSE]

print(table(meta$condition, meta$sample_type))
print(table(meta$condition, meta$study))

if (all(expr >= 0) && all(abs(expr - round(expr)) < .Machine$double.eps^0.5))
  warning("Expression matrix looks integer-like. Confirm it is truly log2(norm_count + 1).")

# ------------------------------------------------------------------
# 6) Expressao diferencial com limma em log2(normalizado + 1)
#    limma requer: genes (linhas) x amostras (colunas)
# ------------------------------------------------------------------
E <- t(expr)
storage.mode(E) <- "numeric"

design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)

cat("dim(expr)   [amostras x genes] =", paste(dim(expr),   collapse = " x "), "\n")
cat("dim(E)      [genes x amostras] =", paste(dim(E),      collapse = " x "), "\n")
cat("dim(design) [amostras x coef]  =", paste(dim(design), collapse = " x "), "\n")
stopifnot(ncol(E) == nrow(design))

fit  <- lmFit(E, design)
fit2 <- contrasts.fit(
  fit,
  makeContrasts(LIHC_vs_Normal = LIHC - Normal, levels = design)
)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2, coef = "LIHC_vs_Normal",
  number = Inf, adjust.method = "BH", sort.by = "P"
) |>
  tibble::rownames_to_column(var = "gene_symbol")

# Classificacao dos DEGs
deg$regulation <- dplyr::case_when(
  deg$adj.P.Val < 0.05 & deg$logFC >  1 ~ "Up_LIHC",
  deg$adj.P.Val < 0.05 & deg$logFC < -1 ~ "Down_LIHC",
  TRUE ~ "NS"
)

print(table(deg$regulation))
print(head(deg, 20))

rio::export(deg, "DEG_LIHC_vs_Normal_limma_log2normcount.csv")

# Genes KEGG ausentes nos resultados de DE
missing_kegg_genes <- setdiff(genes_via, deg$gene_symbol)
cat("KEGG genes nao encontrados nos DEGs:", paste(missing_kegg_genes, collapse = ", "), "\n")

# ------------------------------------------------------------------
# 7) Volcano plot
# ------------------------------------------------------------------
lfc_cutoff <- 1
fdr_cutoff <- 0.05

deg$volcano_class <- dplyr::case_when(
  deg$adj.P.Val < fdr_cutoff & deg$logFC >  lfc_cutoff ~ "Up",
  deg$adj.P.Val < fdr_cutoff & deg$logFC < -lfc_cutoff ~ "Down",
  TRUE ~ "NS"
)

p_volcano <- ggplot(deg, aes(logFC, -log10(adj.P.Val), color = volcano_class)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(Up = "blue", Down = "purple", NS = "grey70")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(
    data = subset(deg, volcano_class %in% c("Up", "Down")),
    aes(label = gene_symbol),
    size = 3.8, max.overlaps = Inf,
    box.padding = 0.4, point.padding = 0.3, segment.color = "grey50"
  ) +
  labs(
    title = "Genes up e downregulados - Via LA em LIHC",
    x     = "log2 Fold Change",
    y     = "-log10(FDR)",
    color = "Regulation"
  ) +
  theme_classic(base_size = 14)

print(p_volcano)
ggsave("outputs/Volcano_LIHC_vs_Normal.png", p_volcano,
       width = 7, height = 6, dpi = 300)

# ------------------------------------------------------------------
# 8) Rede KEGG - DEGs na via de lipideos e aterosclerose
# ------------------------------------------------------------------

# Definir genes Up/Down
up_genes   <- deg |> dplyr::filter(regulation == "Up_LIHC")   |> dplyr::pull(gene_symbol) |> unique()
down_genes <- deg |> dplyr::filter(regulation == "Down_LIHC") |> dplyr::pull(gene_symbol) |> unique()
deg_genes  <- sort(unique(c(up_genes, down_genes)))

cat("Up genes:",   length(up_genes),   "\n")
cat("Down genes:", length(down_genes), "\n")
cat("Total DEGs:", length(deg_genes),  "\n")

# Genes DEG que estão na via KEGG hsa05417
deg_in_pathway <- base::intersect(deg_genes, genes_via)
cat("DEGs na via KEGG hsa05417:", length(deg_in_pathway), "\n")

if (length(deg_in_pathway) < 2) 
  stop("DEGs insuficientes na via KEGG para construir rede.")

# Criar anotações de regulação
node_annot_kegg <- tibble(
  gene_symbol = deg_in_pathway,
  regulation  = dplyr::case_when(
    gene_symbol %in% up_genes   ~ "Up",
    gene_symbol %in% down_genes ~ "Down",
    TRUE ~ "NS"
  )
)

# Criar rede KEGG com genes da via
# Se dois genes estão na mesma via, são conectados
edges_kegg <- tibble()
for (i in 1:(length(deg_in_pathway) - 1)) {
  for (j in (i+1):length(deg_in_pathway)) {
    edges_kegg <- dplyr::bind_rows(
      edges_kegg,
      tibble(
        from   = deg_in_pathway[i],
        to     = deg_in_pathway[j],
        weight = 900  # Conexão pela via KEGG hsa05417
      )
    )
  }
}

if (nrow(edges_kegg) == 0)
  stop("Nenhuma aresta gerada para rede KEGG.")

cat("Arestas KEGG:", nrow(edges_kegg), "\n")

# Construir grafo igraph
g_kegg <- igraph::graph_from_data_frame(
  d        = edges_kegg |> dplyr::transmute(from, to, weight),
  directed = FALSE,
  vertices = node_annot_kegg |> dplyr::transmute(name = gene_symbol, regulation)
)

# Maior componente conexo
comps_kegg <- igraph::components(g_kegg)
giant_kegg <- which.max(comps_kegg$csize)
g_cc_kegg  <- igraph::induced_subgraph(g_kegg, igraph::V(g_kegg)[comps_kegg$membership == giant_kegg])
g_cc_kegg  <- igraph::delete_vertices(g_cc_kegg, igraph::V(g_cc_kegg)[igraph::degree(g_cc_kegg) == 0])

cat("Nos (grafo KEGG):", igraph::vcount(g_kegg),   "| Arestas:", igraph::ecount(g_kegg),   "\n")
cat("Nos (maior CC):   ", igraph::vcount(g_cc_kegg), "| Arestas:", igraph::ecount(g_cc_kegg), "\n")

# Plot Rede KEGG 2D (Kamada-Kawai layout)
set.seed(4721)

p_kegg_net <- ggraph(g_cc_kegg, layout = "kk") +
  geom_edge_link(aes(width = weight), alpha = 0.35) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  geom_node_point(aes(color = regulation), size = 4) +
  scale_color_manual(values = c(Up = "blue", Down = "purple", NS = "grey75")) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
  theme_void(base_size = 14) +
  labs(title = "KEGG Network - Via hsa05417 (Lipideos e Aterosclerose)")

print(p_kegg_net)
ggsave("outputs/Network_KEGG_DEGs_hsa05417.png", p_kegg_net,
       width = 9, height = 7, dpi = 300)

# Plot Rede KEGG 3D (plotly)
lay_3d_kegg <- igraph::layout_with_fr(g_cc_kegg, dim = 3)
igraph::V(g_cc_kegg)$x <- lay_3d_kegg[, 1]
igraph::V(g_cc_kegg)$y <- lay_3d_kegg[, 2]
igraph::V(g_cc_kegg)$z <- lay_3d_kegg[, 3]

edge_width_3d_kegg <- scales::rescale(igraph::E(g_cc_kegg)$weight, to = c(2, 10))
node_colors_kegg   <- ifelse(igraph::V(g_cc_kegg)$regulation == "Up",   "#D55E00",
                      ifelse(igraph::V(g_cc_kegg)$regulation == "Down",  "#0072B2", "grey75"))
node_size_3d_kegg  <- scales::rescale(igraph::degree(g_cc_kegg), to = c(6, 16))

edge_list_kegg <- igraph::as_edgelist(g_cc_kegg)
p_3d_kegg      <- plot_ly()

for (i in seq_len(nrow(edge_list_kegg))) {
  v0   <- edge_list_kegg[i, 1]
  v1   <- edge_list_kegg[i, 2]
  p_3d_kegg <- add_trace(
    p_3d_kegg,
    x = c(igraph::V(g_cc_kegg)[v0]$x, igraph::V(g_cc_kegg)[v1]$x),
    y = c(igraph::V(g_cc_kegg)[v0]$y, igraph::V(g_cc_kegg)[v1]$y),
    z = c(igraph::V(g_cc_kegg)[v0]$z, igraph::V(g_cc_kegg)[v1]$z),
    type = "scatter3d", mode = "lines",
    line = list(color = "#4A4A4A", width = edge_width_3d_kegg[i]),
    hoverinfo = "none", showlegend = FALSE
  )
}

p_3d_kegg <- add_trace(
  p_3d_kegg,
  x    = igraph::V(g_cc_kegg)$x,
  y    = igraph::V(g_cc_kegg)$y,
  z    = igraph::V(g_cc_kegg)$z,
  type = "scatter3d", mode = "markers",
  marker    = list(size = node_size_3d_kegg, color = node_colors_kegg,
                   line = list(color = "black", width = 0.8)),
  text      = igraph::V(g_cc_kegg)$name,
  hoverinfo = "text", showlegend = FALSE
)

p_3d_kegg <- plotly::layout(
  p_3d_kegg,
  scene = list(
    xaxis = list(showticklabels = FALSE, title = ""),
    yaxis = list(showticklabels = FALSE, title = ""),
    zaxis = list(showticklabels = FALSE, title = "")
  ),
  title = "KEGG Network 3D - Via hsa05417 (Lipideos e Aterosclerose)"
)

p_3d_kegg

# ------------------------------------------------------------------
# 9) Mapeamento de genes para ENTREZ ID (para enriquecimento KEGG)
# ------------------------------------------------------------------
# Esta seção prepara os IDs necessários para análise KEGG abaixo
# (Os genes Up e Down já estão mapeados na seção 10)

# ------------------------------------------------------------------
# 10) Enriquecimento KEGG - clusterProfiler
# ------------------------------------------------------------------

# Converter simbolos genicas -> ENTREZ ID
map_to_entrez <- function(symbols) {
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db, drop = TRUE) |>
    dplyr::pull(ENTREZID)
}

entrez_up   <- map_to_entrez(up_genes)
entrez_down <- map_to_entrez(down_genes)

cat("Up genes com ENTREZ ID:",   length(entrez_up),   "\n")
cat("Down genes com ENTREZ ID:", length(entrez_down), "\n")

kegg_up <- enrichKEGG(
  gene         = entrez_up,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_down <- enrichKEGG(
  gene         = entrez_down,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Exportar tabelas KEGG
if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
  rio::export(as.data.frame(kegg_up), "outputs/KEGG_enrichment_Up.csv")
  cat("KEGG Up - vias significativas:", sum(kegg_up@result$p.adjust < 0.05), "\n")
}
if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
  rio::export(as.data.frame(kegg_down), "outputs/KEGG_enrichment_Down.csv")
  cat("KEGG Down - vias significativas:", sum(kegg_down@result$p.adjust < 0.05), "\n")
}

# Dotplot KEGG combinado
kegg_up_df   <- as.data.frame(kegg_up)   |> dplyr::mutate(set_label = "Up_LIHC")
kegg_down_df <- as.data.frame(kegg_down) |> dplyr::mutate(set_label = "Down_LIHC")

kegg_plot_df <- dplyr::bind_rows(kegg_up_df, kegg_down_df) |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::mutate(
    GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x))),
    logFDR        = -log10(p.adjust),
    Description   = stringr::str_wrap(Description, width = 40)
  ) |>
  dplyr::group_by(set_label) |>
  dplyr::slice_max(order_by = logFDR, n = 15) |>
  dplyr::ungroup()

if (nrow(kegg_plot_df) > 0) {

  kegg_plot_df$Description <- factor(
    kegg_plot_df$Description,
    levels = rev(unique(kegg_plot_df$Description[order(kegg_plot_df$logFDR)]))
  )

  p_kegg <- ggplot(kegg_plot_df,
                   aes(x = logFDR, y = Description,
                       size = Count, color = logFDR)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ set_label, scales = "free_y") +
    scale_color_gradient(low = "#74add1", high = "#d73027",
                         name = "-log10(FDR)") +
    scale_size(range = c(3, 9), name = "Gene count") +
    labs(
      x     = "-log10(FDR)",
      y     = "KEGG Pathway",
      title = "KEGG Pathway Enrichment - clusterProfiler"
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text  = element_text(size = 12, face = "bold"),
      plot.title  = element_text(hjust = 0.5)
    )

  print(p_kegg)
  ggsave("outputs/KEGG_dotplot_enrichment.png", p_kegg,
         width = 12, height = 8, dpi = 300)
} else {
  message("Nenhuma via KEGG significativa para plotar.")
}



# ------------------------------------------------------------------
# 13) Heatmap - Top 20 DEGs (ggplot2)
# ------------------------------------------------------------------
top20     <- deg |> dplyr::arrange(adj.P.Val) |> dplyr::slice_head(n = 20)
top_genes <- base::intersect(top20$gene_symbol, rownames(E))

mat_top    <- E[top_genes, ]
mat_scaled <- t(scale(t(mat_top)))   # z-score por gene

ord           <- order(meta$condition)
mat_scaled    <- mat_scaled[, ord]
condition_ord <- meta$condition[ord]

df_long <- as.data.frame(mat_scaled) |>
  dplyr::mutate(Gene = rownames(mat_scaled)) |>
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Zscore")

df_long$Condition <- factor(
  condition_ord[match(df_long$Sample, colnames(mat_scaled))]
)
df_long$Gene <- factor(df_long$Gene, levels = rev(top20$gene_symbol))

p_heatmap <- ggplot(df_long, aes(x = Sample, y = Gene, fill = Zscore)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-2, 2), name = "Z-score"
  ) +
  facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),
    strip.text      = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_heatmap)
ggsave("outputs/Heatmap_Top20_DEGs.png", p_heatmap,
       width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------------
# 14) Sankey plot - Top 20 DEGs por |logFC| (somente Up/Down)
# ------------------------------------------------------------------
sankey_df <- deg |>
  dplyr::filter(regulation != "NS") |>
  dplyr::mutate(
    vol_class = ifelse(logFC > 0, "Up", "Down"),
    abs_logFC = abs(logFC)
  ) |>
  dplyr::arrange(dplyr::desc(abs_logFC)) |>
  dplyr::slice_head(n = 20) |>
  dplyr::mutate(gene_symbol = factor(gene_symbol, levels = gene_symbol))

p_sankey <- ggplot(sankey_df,
                   aes(axis1 = vol_class, axis2 = gene_symbol, y = abs_logFC)) +
  geom_alluvium(aes(fill = vol_class),
                width = 0.28, alpha = 0.65, curve_type = "cubic") +
  geom_stratum(width = 0.32, fill = "white",
               color = "#E8E8E8", linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.3, color = "#2C3E50") +
  scale_fill_manual(values = c(Up = "#C9B2B2", Down = "#9FB7B7")) +
  scale_x_discrete(limits = c("Regulation", "Gene"),
                   expand = c(0.08, 0.05)) +
  ylab("log2 Fold Change") +
  xlab(NULL) +
  theme_minimal() +
  theme(
    legend.position  = "none",
    panel.grid       = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.y      = element_blank(),
    axis.text.x      = element_text(size = 10, color = "#4A4A4A",
                                    margin = margin(t = 5)),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(20, 20, 15, 20)
  )

print(p_sankey)
ggsave("outputs/Sankey_plot.png", p_sankey,
       width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------------
# 15) PCA - Todos os genes da via (presentes no DEG)
# ------------------------------------------------------------------
genes_validos_pca1 <- base::intersect(deg$gene_symbol, rownames(E))
mat_pca1           <- E[genes_validos_pca1, ]
pca_all            <- prcomp(t(mat_pca1), scale. = TRUE)

pca_df1 <- data.frame(
  PC1       = pca_all$x[, 1],
  PC2       = pca_all$x[, 2],
  Condition = meta$condition
)

p_pca1 <- ggplot(pca_df1, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 0.7) +
  labs(
    title = "PCA - Todos os genes da via",
    x = paste0("PC1 (", round(summary(pca_all)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_all)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

print(p_pca1)
ggsave("outputs/PCA_All_Pathway_Genes.png", p_pca1,
       width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------
# 16) PCA - DEGs da via hsa05417
# ------------------------------------------------------------------
genes_via_deg  <- base::intersect(genes_via, deg$gene_symbol)
genes_validos2 <- base::intersect(genes_via_deg, rownames(E))

if (length(genes_validos2) < 2)
  stop("Menos de 2 genes da via encontrados na matriz de expressao.")

mat_pca2 <- E[genes_validos2, , drop = FALSE]
mat_pca2 <- mat_pca2[apply(mat_pca2, 1, var) > 0, , drop = FALSE]

pca_via <- prcomp(t(mat_pca2), scale. = TRUE)

# Criar pca_df2 especifico para esta PCA (nao reutilizar pca_df1)
pca_df2 <- data.frame(
  PC1       = pca_via$x[, 1],
  PC2       = pca_via$x[, 2],
  Condition = meta$condition
)

p_pca2 <- ggplot(pca_df2, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 0.7) +
  labs(
    title = "PCA - DEGs da via hsa05417",
    x = paste0("PC1 (", round(summary(pca_via)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_via)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

print(p_pca2)
ggsave("outputs/PCA_DE_Genes.png", p_pca2,
       width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------
# 17) Heatmap geral - todos os DEGs (pheatmap)
# ------------------------------------------------------------------
all_deg_genes <- base::intersect(
  deg$gene_symbol[deg$regulation != "NS"],
  rownames(E)
)

mat_all        <- E[all_deg_genes, ]
mat_all_scaled <- t(scale(t(mat_all)))

# Clamp outliers para paleta simetrica
mat_all_scaled <- pmin(pmax(mat_all_scaled, -3), 3)

# Ordenar amostras por condicao
ord_all            <- order(meta$condition)
mat_all_scaled_ord <- mat_all_scaled[, ord_all]

# annotation_col: row.names DEVE bater com colnames da matriz
annotation_samples <- data.frame(
  Condition   = meta$condition[ord_all],
  Sample_Type = meta$sample_type[ord_all],
  Study       = meta$study[ord_all],
  row.names   = colnames(mat_all_scaled_ord)
)

pheatmap(
  mat_all_scaled_ord,
  scale             = "none",          # ja escalado acima
  annotation_col    = annotation_samples,
  clustering_method = "ward.D2",
  main              = "Heatmap - Todos os DEGs (z-score)",
  color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
  breaks            = seq(-3, 3, length.out = 51),
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  fontsize_row      = 8,
  fontsize_col      = 8,
  height            = 10,
  width             = 10,
  filename          = "outputs/Heatmap_All_DEGs.png"
)

# ------------------------------------------------------------------
# Salvar workspace
# ------------------------------------------------------------------
save.image(file = "workspace_LIHC.RData")
gc()
