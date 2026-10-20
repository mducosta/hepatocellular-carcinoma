# =============================================================================
# 01_analise_exploratoria.R
# -----------------------------------------------------------------------------
# Análise exploratória:
#   - fenótipos clínicos por grupo (sexo, idade, estágio, grau);
#   - distribuição das amostras;
#   - UMAP (condição, estudo, sexo);
#   - heatmap dos genes mais variáveis;
#   - expressão de marcadores conhecidos de HCC.
#
# Saídas: resultados/tabelas/*.csv e resultados/figuras/*.png
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
})

set.seed(2026)  # reprodutibilidade dos jitters/posições

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_proc <- file.path(ROOT, "dados", "processed")
dir_tab  <- file.path(ROOT, "resultados", "tabelas")
dir_fig  <- file.path(ROOT, "resultados", "figuras")
dir_qc   <- file.path(ROOT, "resultados", "qc")
for (d in c(dir_tab, dir_fig, dir_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 0. Tema e utilidades
# -----------------------------------------------------------------------------
theme_seminario <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, colour = "grey30"),
      axis.text = element_text(colour = "grey20")
    )
}

cores_grupo <- c(
  "TCGA_Tumor"    = "#C0392B",
  "TCGA_Adjacent" = "#E67E22",
  "GTEx_Normal"   = "#2471A3",
  "TCGA_Recurrent"= "#7D3C98"
)
cores_sexo <- c("female" = "#D35400", "male" = "#2E86C1")

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-16) return("< 2.2e-16")
  formatC(p, format = "e", digits = 2)
}

# -----------------------------------------------------------------------------
# 1. Dados
# -----------------------------------------------------------------------------
meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
              sep = "\t", data.table = FALSE)
meta$condition <- factor(meta$condition, levels = c("Tumor", "Adjacent", "Normal"))
meta$group <- factor(meta$group, levels = c("TCGA_Tumor", "TCGA_Adjacent", "GTEx_Normal"))
meta$sex <- factor(meta$sex, levels = c("female", "male"))
meta$stage <- factor(meta$stage, levels = c("I", "II", "III", "IV"))
meta$grade <- factor(meta$grade, levels = c("G1", "G2", "G3", "G4"))

expr_file <- file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz")
if (!file.exists(expr_file)) {
  expr_file <- file.path(dir_proc, "expressao_LA_pathway.tsv.gz")
  warning("Matriz completa ausente; usando painel LA.")
}
expr_dt <- fread(expr_file, sep = "\t", data.table = FALSE)
rownames(expr_dt) <- expr_dt$gene
expr <- as.matrix(expr_dt[, setdiff(names(expr_dt), "gene")])
storage.mode(expr) <- "numeric"
expr <- expr[, intersect(meta$sample, colnames(expr)), drop = FALSE]
meta <- meta[match(colnames(expr), meta$sample), ]
meta$condition <- droplevels(meta$condition)
meta$group <- droplevels(meta$group)
cat(sprintf("[DADOS] %d genes x %d amostras\n", nrow(expr), ncol(expr)))

# -----------------------------------------------------------------------------
# 2. Tabelas de fenótipos
# -----------------------------------------------------------------------------
# 2.1 Amostras por grupo (apenas grupos observados; nunca combinações vazias)
tab_amostras <- meta %>%
  filter(!is.na(group)) %>%
  count(group, name = "n") %>%
  filter(n > 0) %>%
  transmute(Grupo = factor(group, levels = levels(meta$group)), n = n)
write.csv(tab_amostras, file.path(dir_tab, "01_amostras_por_grupo.csv"), row.names = FALSE)

# 2.2 Sexo
tab_sexo <- as.data.frame.matrix(table(meta$group, meta$sex))
tab_sexo$Total <- rowSums(tab_sexo)
tab_sexo <- cbind(Grupo = rownames(tab_sexo), tab_sexo)
write.csv(tab_sexo, file.path(dir_tab, "02_sexo_por_grupo.csv"), row.names = FALSE)

# Qui-quadrado: TCGA tumor vs GTEx normal; TCGA tumor vs TCGA adjacente
chisq_sex <- function(g1, g2) {
  sub <- meta[meta$group %in% c(g1, g2), ]
  t <- table(droplevels(sub$group), sub$sex)
  if (any(dim(t) < 2)) return(NA_real_)
  suppressWarnings(chisq.test(t)$p.value)
}
p_sex_tcgavs_gtex <- chisq_sex("TCGA_Tumor", "GTEx_Normal")
p_sex_tcgavs_adj  <- chisq_sex("TCGA_Tumor", "TCGA_Adjacent")
p_sex_gtex_adj    <- chisq_sex("TCGA_Adjacent", "GTEx_Normal")

# 2.3 Idade
tab_idade <- meta %>%
  filter(!is.na(age)) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    media = mean(age),
    dp = sd(age),
    mediana = median(age),
    q1 = quantile(age, 0.25),
    q3 = quantile(age, 0.75),
    min = min(age),
    max = max(age),
    .groups = "drop"
  )
write.csv(tab_idade, file.path(dir_tab, "03_idade_por_grupo.csv"), row.names = FALSE)

p_idade_kw <- if (length(unique(na.omit(meta$group))) >= 2) {
  kruskal.test(age ~ group, data = meta)$p.value
} else NA_real_
p_idade_tcgavs_gtex <- wilcox.test(age ~ group,
  data = subset(meta, group %in% c("TCGA_Tumor", "GTEx_Normal")))$p.value
p_idade_tcgavs_adj <- wilcox.test(age ~ group,
  data = subset(meta, group %in% c("TCGA_Tumor", "TCGA_Adjacent")))$p.value

# 2.4 Estágio e grau (TCGA tumor)
tcga_tumor <- subset(meta, group == "TCGA_Tumor")
tab_estagio <- as.data.frame(table(Estagio = tcga_tumor$stage, useNA = "ifany"))
names(tab_estagio)[1] <- "Estagio"
tab_estagio$percentual <- round(100 * tab_estagio$Freq / sum(tab_estagio$Freq), 1)
write.csv(tab_estagio, file.path(dir_tab, "04_estagio_TCGA_tumor.csv"), row.names = FALSE)

tab_grau <- as.data.frame(table(Grau = tcga_tumor$grade, useNA = "ifany"))
names(tab_grau)[1] <- "Grau"
tab_grau$percentual <- round(100 * tab_grau$Freq / sum(tab_grau$Freq), 1)
write.csv(tab_grau, file.path(dir_tab, "05_grau_TCGA_tumor.csv"), row.names = FALSE)

tab_estagio_grau <- as.data.frame.matrix(table(tcga_tumor$stage_group, tcga_tumor$grade_group, useNA = "ifany"))
tab_estagio_grau <- cbind(Estagio = rownames(tab_estagio_grau), tab_estagio_grau)
write.csv(tab_estagio_grau, file.path(dir_tab, "06_estagio_x_grau_TCGA.csv"), row.names = FALSE)

# 2.5 Resumo clínico TCGA (tumor)
clin_resumo <- data.frame(
  variavel = c("n", "idade_media_dp", "idade_mediana", "sexo_feminino_n",
               "sexo_masculino_n", "estagio_I", "estagio_II", "estagio_III",
               "estagio_IV", "grau_G1", "grau_G2", "grau_G3", "grau_G4",
               "AFP_mediana", "fibrose_Ishak_mediana"),
  valor = c(
    nrow(tcga_tumor),
    paste0(round(mean(tcga_tumor$age, na.rm = TRUE), 1), " +/- ",
           round(sd(tcga_tumor$age, na.rm = TRUE), 1)),
    round(median(tcga_tumor$age, na.rm = TRUE), 1),
    sum(tcga_tumor$sex == "female", na.rm = TRUE),
    sum(tcga_tumor$sex == "male", na.rm = TRUE),
    sum(tcga_tumor$stage == "I", na.rm = TRUE),
    sum(tcga_tumor$stage == "II", na.rm = TRUE),
    sum(tcga_tumor$stage == "III", na.rm = TRUE),
    sum(tcga_tumor$stage == "IV", na.rm = TRUE),
    sum(tcga_tumor$grade == "G1", na.rm = TRUE),
    sum(tcga_tumor$grade == "G2", na.rm = TRUE),
    sum(tcga_tumor$grade == "G3", na.rm = TRUE),
    sum(tcga_tumor$grade == "G4", na.rm = TRUE),
    round(median(tcga_tumor$afp, na.rm = TRUE), 2),
    round(median(tcga_tumor$fibrosis_ishak, na.rm = TRUE), 2)
  )
)
write.csv(clin_resumo, file.path(dir_tab, "07_resumo_clinico_TCGA.csv"), row.names = FALSE)

# 2.6 Testes de fenótipos
testes <- data.frame(
  comparacao = c("Sexo: TCGA tumor vs GTEx normal",
                 "Sexo: TCGA tumor vs TCGA adjacente",
                 "Sexo: TCGA adjacente vs GTEx normal",
                 "Idade: TCGA tumor vs GTEx normal (Wilcoxon)",
                 "Idade: TCGA tumor vs TCGA adjacente (Wilcoxon)",
                 "Idade: Kruskal-Wallis entre os 3 grupos"),
  p_valor = c(p_sex_tcgavs_gtex, p_sex_tcgavs_adj, p_sex_gtex_adj,
              p_idade_tcgavs_gtex, p_idade_tcgavs_adj, p_idade_kw)
)
testes$p_ajustado_BH <- p.adjust(testes$p_valor, method = "BH")
testes$significativo <- ifelse(testes$p_ajustado_BH < 0.05, "Sim", "Não")
write.csv(testes, file.path(dir_tab, "08_testes_fenotipos.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Figuras de fenótipos
# -----------------------------------------------------------------------------
# 3.1 Amostras (barras apenas para grupos com n > 0)
p1 <- ggplot(tab_amostras, aes(x = Grupo, y = n, fill = Grupo)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.35, size = 4) +
  scale_fill_manual(values = cores_grupo) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Amostras por grupo (análise principal)",
       subtitle = paste0("TCGA-LIHC (tumor primário e adjacente) e GTEx (fígado normal); ",
                         "apenas grupos com n > 0"),
       x = NULL, y = "Número de amostras", fill = "Grupo") +
  theme_seminario() +
  theme(legend.position = "none")
ggsave(file.path(dir_fig, "01_amostras_por_grupo.png"), p1,
       width = 7, height = 4.5, dpi = 300)

# 3.2 Sexo
p2 <- meta %>%
  filter(!is.na(sex)) %>%
  count(group, sex) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = group, y = prop, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = cores_sexo) +
  labs(title = "Composição por sexo",
       subtitle = paste0("Qui-quadrado TCGA tumor vs GTEx normal: p = ",
                         fmt_p(p_sex_tcgavs_gtex)),
       x = NULL, y = "Proporção", fill = "Sexo") +
  theme_seminario()
ggsave(file.path(dir_fig, "02_sexo_por_grupo.png"), p2,
       width = 7, height = 4.5, dpi = 300)

# 3.3 Idade
p3 <- meta %>%
  filter(!is.na(age)) %>%
  ggplot(aes(x = group, y = age, fill = group)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9, colour = "grey20") +
  geom_jitter(width = 0.10, size = 0.9, alpha = 0.35, colour = "grey20") +
  scale_fill_manual(values = cores_grupo) +
  labs(title = "Distribuição de idade por grupo (violino + boxplot)",
       subtitle = paste0("Wilcoxon TCGA tumor vs GTEx normal: p = ",
                         fmt_p(p_idade_tcgavs_gtex)),
       x = NULL, y = "Idade (anos)", fill = "Grupo") +
  theme_seminario() +
  theme(legend.position = "none")
ggsave(file.path(dir_fig, "03_idade_por_grupo.png"), p3,
       width = 6.5, height = 4.5, dpi = 300)

# 3.4 Estágio
p4 <- tcga_tumor %>%
  filter(!is.na(stage)) %>%
  count(stage) %>%
  ggplot(aes(x = stage, y = n, fill = stage)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 4) +
  scale_fill_brewer(palette = "OrRd") +
  labs(title = "Estágio AJCC — TCGA-LIHC (tumor primário)",
       subtitle = paste0("n = ", sum(!is.na(tcga_tumor$stage)),
                         " tumores com estágio informado"),
       x = "Estágio", y = "Número de tumores", fill = "Estágio") +
  theme_seminario() +
  theme(legend.position = "none")
ggsave(file.path(dir_fig, "04_estagio_TCGA.png"), p4,
       width = 6, height = 4.2, dpi = 300)

# 3.5 Grau
p5 <- tcga_tumor %>%
  filter(!is.na(grade)) %>%
  count(grade) %>%
  ggplot(aes(x = grade, y = n, fill = grade)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 4) +
  scale_fill_brewer(palette = "Purples") +
  labs(title = "Grau histológico — TCGA-LIHC (tumor primário)",
       subtitle = paste0("n = ", sum(!is.na(tcga_tumor$grade)),
                         " tumores com grau informado"),
       x = "Grau", y = "Número de tumores", fill = "Grau") +
  theme_seminario() +
  theme(legend.position = "none")
ggsave(file.path(dir_fig, "05_grau_TCGA.png"), p5,
       width = 6, height = 4.2, dpi = 300)

# 3.6 Estágio x grau
tab_eg <- as.data.frame(table(Estagio = tcga_tumor$stage_group,
                              Grau = tcga_tumor$grade_group))
p5b <- ggplot(tab_eg, aes(x = Estagio, y = Freq, fill = Grau)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("G1-G2" = "#9BC2E6", "G3-G4" = "#1F4E79")) +
  labs(title = "Estágio x grau histológico — TCGA-LIHC",
       x = "Estágio agrupado", y = "Número de tumores", fill = "Grau") +
  theme_seminario()
ggsave(file.path(dir_fig, "06_estagio_x_grau_TCGA.png"), p5b,
       width = 6.5, height = 4.2, dpi = 300)

# -----------------------------------------------------------------------------
# 4. UMAP (substitui a PCA)
# -----------------------------------------------------------------------------
# Top 2.000 genes por variância
vars <- apply(expr, 1, var, na.rm = TRUE)
top_var <- names(sort(vars, decreasing = TRUE))[seq_len(min(2000, length(vars)))]
X <- t(expr[top_var, , drop = FALSE])
X <- X[complete.cases(X), , drop = FALSE]

set.seed(2026)
umap_res <- uwot::umap(
  X,
  n_neighbors = 30,
  min_dist = 0.25,
  metric = "cosine",
  scale = TRUE,
  pca = 50,
  seed = 2026,
  n_threads = 1,
  n_sgd_threads = 1,
  verbose = FALSE
)
umap_df <- data.frame(
  sample = rownames(X),
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2]
)
umap_df <- merge(umap_df, meta, by = "sample", sort = FALSE)
write.csv(umap_df, file.path(dir_tab, "09_UMAP_coordenadas.csv"), row.names = FALSE)
write.csv(
  data.frame(
    parametro = c("n_genes", "n_neighbors", "min_dist", "metric", "scale", "pca", "seed"),
    valor = c(length(top_var), 30, 0.25, "cosine", TRUE, 50, 2026)
  ),
  file.path(dir_tab, "10_UMAP_parametros.csv"), row.names = FALSE
)

umap_plot <- function(df, colour, title, subtitle, pal) {
  ggplot(df, aes(UMAP1, UMAP2, colour = .data[[colour]])) +
    geom_point(size = 1.8, alpha = 0.75) +
    stat_ellipse(aes(group = .data[[colour]]), type = "norm",
                 level = 0.68, linewidth = 0.4, show.legend = FALSE) +
    scale_colour_manual(values = pal) +
    labs(title = title, subtitle = subtitle,
         x = "UMAP 1", y = "UMAP 2", colour = colour) +
    theme_seminario()
}

p6 <- umap_plot(umap_df, "group",
               "UMAP das amostras hepáticas",
               paste0("Top 2.000 genes mais variáveis | ", nrow(umap_df), " amostras"),
               cores_grupo)
ggsave(file.path(dir_fig, "07_UMAP_grupo.png"), p6, width = 7, height = 5.5, dpi = 300)

p7 <- umap_plot(umap_df, "study",
               "UMAP por estudo (TCGA vs GTEx)",
               "Separação global indica forte efeito de lote/estudo",
               c("TCGA" = "#C0392B", "GTEX" = "#2471A3"))
ggsave(file.path(dir_fig, "08_UMAP_estudo.png"), p7, width = 7, height = 5.5, dpi = 300)

p8 <- umap_plot(umap_df[umap_df$group == "TCGA_Tumor", ], "sex",
               "UMAP dentro de TCGA-LIHC (tumor)",
               "Coloração por sexo",
               cores_sexo)
ggsave(file.path(dir_fig, "09_UMAP_sexo_TCGA.png"), p8, width = 7, height = 5.5, dpi = 300)

# UMAP com estágio e grau (TCGA tumor)
p9 <- ggplot(umap_df[umap_df$group == "TCGA_Tumor" & !is.na(umap_df$stage), ],
             aes(UMAP1, UMAP2, colour = stage, shape = grade)) +
  geom_point(size = 2.2, alpha = 0.8) +
  scale_colour_brewer(palette = "OrRd") +
  labs(title = "UMAP em TCGA-LIHC: estágio e grau",
       subtitle = "Tumores primários",
       x = "UMAP 1", y = "UMAP 2",
       colour = "Estágio", shape = "Grau") +
  theme_seminario()
ggsave(file.path(dir_fig, "10_UMAP_estagio_grau_TCGA.png"), p9,
       width = 7.5, height = 5.5, dpi = 300)

# Associação das coordenadas UMAP com variáveis clínicas (TCGA tumor)
tcga_umap <- umap_df[umap_df$group == "TCGA_Tumor", ]
tcga_umap$stage_num <- as.numeric(tcga_umap$stage)
tcga_umap$grade_num <- as.numeric(tcga_umap$grade)
assoc <- data.frame(
  variavel = c("idade", "estagio", "grau"),
  UMAP1_rho = c(cor(tcga_umap$UMAP1, tcga_umap$age, method = "spearman", use = "complete.obs"),
                cor(tcga_umap$UMAP1, tcga_umap$stage_num, method = "spearman", use = "complete.obs"),
                cor(tcga_umap$UMAP1, tcga_umap$grade_num, method = "spearman", use = "complete.obs")),
  UMAP2_rho = c(cor(tcga_umap$UMAP2, tcga_umap$age, method = "spearman", use = "complete.obs"),
                cor(tcga_umap$UMAP2, tcga_umap$stage_num, method = "spearman", use = "complete.obs"),
                cor(tcga_umap$UMAP2, tcga_umap$grade_num, method = "spearman", use = "complete.obs"))
)
assoc$UMAP1_p <- c(cor.test(tcga_umap$UMAP1, tcga_umap$age, method = "spearman")$p.value,
                   cor.test(tcga_umap$UMAP1, tcga_umap$stage_num, method = "spearman")$p.value,
                   cor.test(tcga_umap$UMAP1, tcga_umap$grade_num, method = "spearman")$p.value)
assoc$UMAP2_p <- c(cor.test(tcga_umap$UMAP2, tcga_umap$age, method = "spearman")$p.value,
                   cor.test(tcga_umap$UMAP2, tcga_umap$stage_num, method = "spearman")$p.value,
                   cor.test(tcga_umap$UMAP2, tcga_umap$grade_num, method = "spearman")$p.value)
write.csv(assoc, file.path(dir_tab, "11_UMAP_associacao_clinica.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Heatmap dos genes mais variáveis
# -----------------------------------------------------------------------------
top50 <- names(sort(vars, decreasing = TRUE))[seq_len(min(50, length(vars)))]
hm <- expr[top50, , drop = FALSE]
hm_z <- t(scale(t(hm)))
hm_z[hm_z > 2] <- 2
hm_z[hm_z < -2] <- -2

annot_col <- data.frame(
  Grupo = meta$group,
  Sexo = meta$sex,
  row.names = meta$sample
)
annot_colors <- list(
  Grupo = cores_grupo,
  Sexo = cores_sexo
)
png(file.path(dir_fig, "11_heatmap_top50_variaveis.png"),
    width = 11, height = 9, units = "in", res = 300)
pheatmap(hm_z,
         annotation_col = annot_col,
         annotation_colors = annot_colors,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         fontsize_row = 7,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Top 50 genes mais variáveis (z-score)")
dev.off()

# -----------------------------------------------------------------------------
# 6. Marcadores conhecidos de HCC
# -----------------------------------------------------------------------------
marcadores <- c("AFP", "GPC3", "SPINK1", "BIRC5", "EPCAM", "KRT19",
                "VIM", "CDH1", "MMP9", "COL1A1", "SPP1", "CXCL8",
                "TGFB1", "VEGFA", "MYC", "CTNNB1", "GLUL", "LGR5")
marcadores <- intersect(marcadores, rownames(expr))
if (length(marcadores) > 0) {
  mar_df <- as.data.frame(t(expr[marcadores, , drop = FALSE]))
  mar_df$sample <- rownames(mar_df)
  mar_long <- tidyr::pivot_longer(mar_df, cols = all_of(marcadores),
                                  names_to = "gene", values_to = "expressao")
  mar_long <- merge(mar_long, meta[, c("sample", "group")], by = "sample")
  mar_long$gene <- factor(mar_long$gene, levels = marcadores)
  p_mar <- ggplot(mar_long, aes(x = group, y = expressao, fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.55, colour = NA) +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9, colour = "grey20") +
    facet_wrap(~ gene, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = cores_grupo) +
    labs(title = "Expressão de marcadores de HCC e remodelamento tecidual",
         subtitle = "log2(TPM + 0,001); violino + boxplot",
         x = NULL, y = "Expressão (log2 TPM)", fill = "Grupo") +
    theme_seminario(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "bottom")
  ggsave(file.path(dir_fig, "12_violino_marcadores.png"), p_mar,
         width = 13, height = 8, dpi = 300)
  write.csv(mar_long, file.path(dir_tab, "12_expressao_marcadores.csv"), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 7. Salvar resumo
# -----------------------------------------------------------------------------
resumo <- c(
  "ANÁLISE EXPLORATÓRIA — RESUMO",
  "=============================",
  paste("Genes analisados:", nrow(expr)),
  paste("Amostras:", ncol(expr)),
  "UMAP: top 2.000 genes variáveis, n_neighbors = 30, min_dist = 0,25, métrica = cosine",
  paste("Genes usados no UMAP:", length(top_var), "| Amostras:", nrow(umap_df)),
  "",
  "Amostras por grupo:",
  paste(capture.output(print(table(meta$group))), collapse = "\n"),
  "",
  "Sexo por grupo:",
  paste(capture.output(print(table(meta$group, meta$sex))), collapse = "\n"),
  "",
  "Testes de fenótipos:",
  paste(capture.output(print(testes, row.names = FALSE)), collapse = "\n"),
  "",
  "Estágio (TCGA tumor):",
  paste(capture.output(print(table(tcga_tumor$stage))), collapse = "\n"),
  "Grau (TCGA tumor):",
  paste(capture.output(print(table(tcga_tumor$grade))), collapse = "\n"),
  "",
  "Associação das coordenadas UMAP com clínica (TCGA tumor):",
  paste(capture.output(print(assoc, row.names = FALSE)), collapse = "\n"),
  "",
  "Sessão R:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(resumo, file.path(dir_qc, "resumo_analise_exploratoria.txt"))
cat(paste(resumo, collapse = "\n"), "\n")
