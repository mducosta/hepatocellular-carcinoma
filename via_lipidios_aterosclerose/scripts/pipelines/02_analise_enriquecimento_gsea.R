# ===================================================================
# ANÁLISE DE ENRIQUECIMENTO FUNCIONAL + GSEA
# LIHC vs Normal — Via Lipid and Atherosclerosis (KEGG hsa05417)
# ===================================================================
# Utiliza como entrada a tabela completa de DEGs já calculada pelo
# pipeline (outputs/deg/DEG_LIHC_vs_Normal_full.csv) e gera:
#   1) ORA (over-representation) com background do genoma completo
#   2) GSEA (rank-based) com genes ordenados por logFC
# Saídas em outputs/enrichment/
# ===================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(rio)
})

PROJECT_ROOT <- normalizePath(".", mustWork = TRUE)
setwd(PROJECT_ROOT)

DEG_FILE <- "outputs/deg/DEG_LIHC_vs_Normal_full.csv"
stopifnot(file.exists(DEG_FILE))

deg <- read.csv(DEG_FILE, stringsAsFactors = FALSE)
cat(sprintf("Genes lidos: %d\n", nrow(deg)))

# ------------------------------------------------------------------
# 0) Mapeamento SYMBOL -> ENTREZID
# ------------------------------------------------------------------
map_sym2entrez <- function(symbols) {
  tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE)
  }, error = function(e) data.frame(SYMBOL = character(), ENTREZID = character()))
}

all_map <- map_sym2entrez(deg$gene_symbol)
cat(sprintf("Genes mapeados para ENTREZ: %d / %d\n", nrow(all_map), nrow(deg)))

deg <- deg %>%
  left_join(all_map, by = c("gene_symbol" = "SYMBOL"))

up_genes   <- deg %>% filter(regulation == "Up_LIHC")   %>% pull(gene_symbol)
down_genes <- deg %>% filter(regulation == "Down_LIHC") %>% pull(gene_symbol)

up_entrez   <- all_map %>% filter(SYMBOL %in% up_genes)   %>% pull(ENTREZID)
down_entrez <- all_map %>% filter(SYMBOL %in% down_genes) %>% pull(ENTREZID)

cat(sprintf("Up: %d genes (%d mapeados) | Down: %d genes (%d mapeados)\n",
            length(up_genes), length(up_entrez),
            length(down_genes), length(down_entrez)))

# ------------------------------------------------------------------
# 1) ORA — enriquecimento com BACKGROUND DO GENOMA COMPLETO
#    (corrige o resultado vazio da versão anterior, que restringia o
#     universo aos 212 genes da via)
# ------------------------------------------------------------------
run_ora <- function(entrez, label) {
  out <- list()
  if (length(entrez) >= 3) {
    out$go <- tryCatch(
      enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) NULL)
    out$kegg <- tryCatch(
      enrichKEGG(gene = entrez, organism = "hsa",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2),
      error = function(e) NULL)
  }
  out
}

ora_up   <- run_ora(up_entrez, "up")
ora_down <- run_ora(down_entrez, "down")

export_ora <- function(obj, path) {
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    rio::export(as.data.frame(obj), path)
    cat(sprintf("  -> %s (%d termos)\n", path, nrow(as.data.frame(obj))))
  } else {
    rio::export(data.frame(), path)
    cat(sprintf("  -> %s (vazio)\n", path))
  }
}

cat("\nORA GO_BP_up:\n"); export_ora(ora_up$go,  "outputs/enrichment/GO_BP_up.csv")
cat("ORA GO_BP_down:\n"); export_ora(ora_down$go,"outputs/enrichment/GO_BP_down.csv")
cat("ORA KEGG_up:\n");   export_ora(ora_up$kegg, "outputs/enrichment/KEGG_up.csv")
cat("ORA KEGG_down:\n"); export_ora(ora_down$kegg,"outputs/enrichment/KEGG_down.csv")

# ------------------------------------------------------------------
# 2) GSEA — genes ordenados por logFC (rank-based)
# ------------------------------------------------------------------
# Lista ranqueada por logFC (todos os 212 genes)
ranked <- deg %>%
  filter(!is.na(ENTREZID)) %>%
  arrange(desc(logFC)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gsea_genes <- ranked$logFC
names(gsea_genes) <- ranked$ENTREZID
cat(sprintf("\nGSEA: lista ranqueada com %d genes\n", length(gsea_genes)))

run_gsea <- function(geneList, fun, ...) {
  tryCatch(fun(geneList = geneList, ...), error = function(e) NULL)
}

# 2.1 GSEA GO (BP)
gsea_go <- run_gsea(gsea_genes, gseGO, ont = "BP",
                    OrgDb = org.Hs.eg.db, eps = 0,
                    pvalueCutoff = 0.1, seed = 4721)
if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  rio::export(as.data.frame(gsea_go), "outputs/enrichment/GSEA_GO_BP.csv")
  cat(sprintf("GSEA GO_BP: %d conjuntos enriquecidos\n", nrow(as.data.frame(gsea_go))))
} else {
  rio::export(data.frame(), "outputs/enrichment/GSEA_GO_BP.csv")
  cat("GSEA GO_BP: vazio\n")
}

# 2.2 GSEA KEGG
gsea_kegg <- run_gsea(gsea_genes, gseKEGG, organism = "hsa",
                      eps = 0, pvalueCutoff = 0.1, seed = 4721)
if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  rio::export(as.data.frame(gsea_kegg), "outputs/enrichment/GSEA_KEGG.csv")
  cat(sprintf("GSEA KEGG: %d conjuntos enriquecidos\n", nrow(as.data.frame(gsea_kegg))))
} else {
  rio::export(data.frame(), "outputs/enrichment/GSEA_KEGG.csv")
  cat("GSEA KEGG: vazio\n")
}

# 2.3 GSEA Hallmark (MSigDB) via fgsea
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

hallmark_list <- split(hallmark$entrez_gene, hallmark$gs_name)
hallmark_list <- lapply(hallmark_list, as.character)

set.seed(4721)
fgsea_res <- fgsea(pathways = hallmark_list, stats = gsea_genes,
                   minSize = 5, maxSize = 500, nPermSimple = 10000)

if (nrow(fgsea_res) > 0) {
  fgsea_out <- fgsea_res %>%
    mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))) %>%
    arrange(padj)
  rio::export(fgsea_out, "outputs/enrichment/GSEA_HALLMARK.csv")
  cat(sprintf("GSEA Hallmark: %d conjuntos testados | %d significativos (padj<0.05)\n",
              nrow(fgsea_out), sum(fgsea_out$padj < 0.05, na.rm = TRUE)))
} else {
  rio::export(data.frame(), "outputs/enrichment/GSEA_HALLMARK.csv")
  cat("GSEA Hallmark: vazio\n")
}

# ------------------------------------------------------------------
# 3) FIGURAS
# ------------------------------------------------------------------

# 3.1 Dotplot ORA (GO BP up/down + KEGG up/down combinados)
ora_plots <- list()
for (nm in c("up", "down")) {
  for (db in c("go", "kegg")) {
    obj <- get(paste0("ora_", nm))[[db]]
    if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
      d <- as.data.frame(obj) %>% head(10) %>% mutate(set = paste0(nm, "_", db))
      ora_plots[[paste0(nm, "_", db)]] <- d
    }
  }
}

if (length(ora_plots) > 0) {
  ora_df <- bind_rows(ora_plots)
  if (nrow(ora_df) > 0) {
    ora_df <- ora_df %>%
      mutate(logFDR = -log10(p.adjust))
    p_ora <- ggplot(ora_df, aes(x = logFDR, y = reorder(Description, logFDR),
                                size = Count, color = logFDR)) +
      geom_point(alpha = 0.9) +
      facet_wrap(~ set, scales = "free_y", ncol = 1) +
      scale_color_gradient(low = "#74add1", high = "#d73027", name = "-log10(FDR)") +
      scale_size(range = c(3, 9), name = "Gene Count") +
      labs(x = "-log10(FDR)", y = "Termo enriquecido",
           title = "Enriquecimento funcional (ORA) — GO BP e KEGG") +
      theme_bw(base_size = 12) +
      theme(axis.text.y = element_text(size = 8),
            strip.text = element_text(size = 10, face = "bold"))
    ggsave("outputs/enrichment/enrichment_dotplot.png", p_ora,
           width = 11, height = 9, dpi = 300)
    cat("\nFigura: outputs/enrichment/enrichment_dotplot.png\n")
  }
}

# 3.2 Barplot GSEA Hallmark (top 15 por |NES|)
if (nrow(fgsea_res) > 0) {
  top_h <- fgsea_res %>%
    arrange(padj) %>%
    filter(!is.na(NES)) %>%
    head(15) %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway),
           pathway = gsub("_", " ", pathway),
           direction = ifelse(NES > 0, "Up_LIHC", "Down_LIHC"))
  
  if (nrow(top_h) > 0) {
    p_h <- ggplot(top_h, aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c(Up_LIHC = "#D73027", Down_LIHC = "#4575B4")) +
      labs(x = NULL, y = "Normalized Enrichment Score (NES)",
           title = "GSEA — Hallmark (MSigDB)",
           subtitle = "Genes da via hsa05417 ordenados por logFC",
           fill = "Direção") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 8))
    ggsave("outputs/enrichment/GSEA_hallmark_barplot.png", p_h,
           width = 10, height = 6, dpi = 300)
    cat("Figura: outputs/enrichment/GSEA_hallmark_barplot.png\n")
  }
}

cat("\n=== ANÁLISE DE ENRIQUECIMENTO + GSEA CONCLUÍDA ===\n")
