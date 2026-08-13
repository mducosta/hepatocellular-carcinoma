# ===================================================================
# ANÁLISE EM 3 GRUPOS — Normal (GTEx) × Adjacente (TCGA) × LIHC
# ===================================================================
# Separa o tecido normal (GTEx), o tecido normal adjacente (TCGA Solid
# Tissue Normal) e o tumor (TCGA Primary Tumor), faz:
#   * QC (PCA por grupo e por estudo)
#   * DE com limma para 3 contrastes
#   * Correção de lote ComBat (opcional) para LIHC × Normal
#   * Análise de sobrevivência (Kaplan-Meier) dos genes-chave
# Saídas em results/3grupos/
# ===================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(limma); library(sva); library(ggplot2); library(rio)
  library(survival); library(survminer)
})

PROJECT_ROOT <- normalizePath(".", mustWork = TRUE)
setwd(PROJECT_ROOT)
dir.create("results/3grupos", showWarnings = FALSE, recursive = TRUE)

DATA_FILE <- "liver.tsv"
stopifnot(file.exists(DATA_FILE))

liver <- read.delim(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
names(liver) <- sub("^_+", "", names(liver))
if ("samples" %in% colnames(liver)) liver$samples <- NULL

# ------------------------------------------------------------------
# Metadados e colunas de gene
# ------------------------------------------------------------------
meta_cols <- c("sample", "detailed_category", "DFI", "DFI.time", "DSS",
               "DSS.time", "gender", "TCGA_GTEX_main_category", "OS",
               "OS.time", "PFI", "PFI.time", "primary disease or tissue",
               "primary_site", "sample_type", "study")

gene_cols <- setdiff(colnames(liver), meta_cols)
gene_cols <- gene_cols[!grepl("^(OS|DSS|DFI|PFI|RFS|age|stage|grade|vital)",
                              gene_cols, ignore.case = TRUE)]

# ------------------------------------------------------------------
# Definir 3 grupos (exclui Recurrent Tumor)
# ------------------------------------------------------------------
liver$group <- dplyr::case_when(
  liver$sample_type == "Normal Tissue"        ~ "Normal",
  liver$sample_type == "Solid Tissue Normal"  ~ "Adjacent",
  liver$sample_type == "Primary Tumor"        ~ "LIHC",
  TRUE ~ "Excluir"
)
liver <- liver[liver$group != "Excluir", ]
liver$group <- factor(liver$group, levels = c("Normal", "Adjacent", "LIHC"))

cat("Grupos:\n")
print(table(liver$group, liver$study))

meta <- liver[, meta_cols, drop = FALSE]
meta$group <- liver$group

# Matriz de expressão (genes x amostras)
E <- t(as.matrix(liver[, gene_cols, drop = FALSE]))
storage.mode(E) <- "numeric"
colnames(E) <- liver$sample
E <- E[apply(E, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
cat(sprintf("\nMatriz: %d genes x %d amostras\n", nrow(E), ncol(E)))

# ------------------------------------------------------------------
# QC — PCA
# ------------------------------------------------------------------
pca <- prcomp(t(E), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Group = meta$group, Study = meta$study)
pv1 <- round(summary(pca)$importance[2, 1] * 100, 1)
pv2 <- round(summary(pca)$importance[2, 2] * 100, 1)

p1 <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c(Normal = "#2E86AB", Adjacent = "#E6A817",
                                LIHC = "#A23B72")) +
  labs(title = "PCA — 3 grupos", subtitle = sprintf("PC1 %s%% | PC2 %s%%", pv1, pv2)) +
  theme_classic(base_size = 14)
ggsave("results/3grupos/PCA_3grupos.png", p1, width = 8, height = 6, dpi = 300)

p2 <- ggplot(pca_df, aes(PC1, PC2, color = Study, shape = Group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  labs(title = "PCA — colorido por estudo (lote)") +
  theme_classic(base_size = 14)
ggsave("results/3grupos/PCA_3grupos_batch.png", p2, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------
# DE com limma — 3 contrastes
# ------------------------------------------------------------------
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)
fit <- lmFit(E, design)

cm <- makeContrasts(
  LIHC_vs_Normal   = LIHC - Normal,
  LIHC_vs_Adjacent = LIHC - Adjacent,
  Adjacent_vs_Normal = Adjacent - Normal,
  levels = design
)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2, robust = TRUE)

results_summary <- list()
for (coef in colnames(cm)) {
  res <- topTable(fit2, coef = coef, number = Inf, sort.by = "P",
                  adjust.method = "BH") |>
    tibble::rownames_to_column(var = "gene_symbol")
  res$regulation <- dplyr::case_when(
    res$adj.P.Val < 0.05 & res$logFC >  1 ~ "Up",
    res$adj.P.Val < 0.05 & res$logFC < -1 ~ "Down",
    TRUE ~ "NS"
  )
  fname <- paste0("results/3grupos/DEG_", coef, ".csv")
  rio::export(res, fname)
  n_up <- sum(res$regulation == "Up")
  n_down <- sum(res$regulation == "Down")
  cat(sprintf("\n%s: %d DEGs (Up=%d Down=%d)\n", coef, n_up + n_down, n_up, n_down))
  results_summary[[coef]] <- c(Up = n_up, Down = n_down,
                               Total = n_up + n_down)
}
rio::export(as.data.frame(do.call(rbind, results_summary)),
            "results/3grupos/resumo_DEGs.csv")

# ------------------------------------------------------------------
# ComBat — correção de lote para LIHC × Normal (TCGA × GTEx)
# ------------------------------------------------------------------
cat("\n--- ComBat (LIHC × Normal, batch = study) ---\n")
keep_ln <- meta$group %in% c("LIHC", "Normal")
meta_ln <- meta[keep_ln, ]
meta_ln$group <- droplevels(meta_ln$group)
E_ln <- E[, meta_ln$sample, drop = FALSE]

E_combat <- tryCatch(
  sva::ComBat(dat = E_ln, batch = meta_ln$study,
              mod = model.matrix(~ meta_ln$group)),
  error = function(e) NULL
)

if (!is.null(E_combat)) {
  design_ln <- model.matrix(~ 0 + group, data = meta_ln)
  colnames(design_ln) <- levels(meta_ln$group)[levels(meta_ln$group) %in% meta_ln$group]
  fit_c <- lmFit(E_combat, design_ln)
  cm_c <- makeContrasts(LIHC_vs_Normal = LIHC - Normal, levels = design_ln)
  fit_c2 <- contrasts.fit(fit_c, cm_c)
  fit_c2 <- eBayes(fit_c2, robust = TRUE)
  res_c <- topTable(fit_c2, coef = 1, number = Inf, sort.by = "P",
                    adjust.method = "BH") |>
    tibble::rownames_to_column(var = "gene_symbol")
  res_c$regulation <- dplyr::case_when(
    res_c$adj.P.Val < 0.05 & res_c$logFC >  1 ~ "Up",
    res_c$adj.P.Val < 0.05 & res_c$logFC < -1 ~ "Down",
    TRUE ~ "NS"
  )
  rio::export(res_c, "results/3grupos/DEG_LIHC_vs_Normal_ComBat.csv")
  cat(sprintf("ComBat: %d DEGs (Up=%d Down=%d)\n",
              sum(res_c$regulation != "NS"),
              sum(res_c$regulation == "Up"),
              sum(res_c$regulation == "Down")))
} else {
  cat("ComBat não executado (erro)\n")
}

# ------------------------------------------------------------------
# Sobrevivência — Kaplan-Meier (apenas LIHC)
# ------------------------------------------------------------------
cat("\n--- Kaplan-Meier (LIHC, OS) ---\n")
lihc_meta <- meta[meta$group == "LIHC", ]
lihc_meta$OS.time <- suppressWarnings(as.numeric(lihc_meta$OS.time))
lihc_meta$OS <- suppressWarnings(as.numeric(lihc_meta$OS))
lihc_meta$PFI.time <- suppressWarnings(as.numeric(lihc_meta$PFI.time))
lihc_meta$PFI <- suppressWarnings(as.numeric(lihc_meta$PFI))

genes_km <- c("HSP90AB1", "MMP9", "CD36", "BAX", "CALML5", "TLR2", "STAT3",
              "MMP1", "ABCA1", "CXCL2")
genes_km <- base::intersect(genes_km, rownames(E))

km_results <- list()
for (g in genes_km) {
  expr_g <- as.numeric(E[g, lihc_meta$sample])
  med <- median(expr_g, na.rm = TRUE)
  grp <- ifelse(expr_g >= med, "Alto", "Baixo")
  df <- data.frame(time = lihc_meta$OS.time, event = lihc_meta$OS, grp = grp)
  df <- df[!is.na(df$time) & df$time > 0, ]
  if (length(unique(df$grp)) < 2 || nrow(df) < 30) next
  fit_km <- survfit(Surv(time, event) ~ grp, data = df)
  lr <- survdiff(Surv(time, event) ~ grp, data = df)
  pval <- 1 - pchisq(lr$chisq, df = 1)
  km_results[[g]] <- data.frame(gene = g, p_logrank = pval, stringsAsFactors = FALSE)
  
  p_km <- ggsurvplot(fit_km, data = df, pval = TRUE, risk.table = TRUE,
                     palette = c("#A23B72", "#2E86AB"),
                     title = paste(g, "— Sobrevida global (LIHC)"))
  ggsave(paste0("results/3grupos/KM_OS_", g, ".png"),
         p_km$plot, width = 7, height = 5, dpi = 300)
}
km_df <- do.call(rbind, km_results)
km_df <- km_df %>% arrange(p_logrank)
rio::export(km_df, "results/3grupos/survival_logrank.csv")
print(km_df)

cat("\n=== ANÁLISE 3 GRUPOS CONCLUÍDA ===\n")
