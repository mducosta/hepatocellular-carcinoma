#!/usr/bin/env Rscript
# =============================================================================
# etapa_29_validacao_externa.R
# Validação externa do sinal tumor vs adjacente (GSE14520) e validação
# interna do núcleo de grau (TCGA, validação cruzada e leave-one-center-out).
# GSE14520 não possui grau de Edmondson; isso é registrado no relatório.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(limma); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))
dir_val <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "validacao_externa")
dir_de  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "expressao_diferencial_covariaveis")
dir.create(dir_val, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
FDR <- 0.05; LFC <- 1

auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]; label <- as.logical(label[ok])
  n1 <- sum(label); n0 <- sum(!label)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
deg_set <- function(tt) tt$gene[tt$regulation != "NS"]

# ---------------------------------------------------------------------------
# GSE14520: tumor vs adjacente
# ---------------------------------------------------------------------------
Xv_dt <- fread(file.path(dir_val, "GSE14520_expression.tsv"), data.table = FALSE)
Xv <- as.matrix(Xv_dt[, -1, drop = FALSE]); rownames(Xv) <- Xv_dt[[1]]; storage.mode(Xv) <- "numeric"
mv <- fread(file.path(dir_val, "GSE14520_metadata.tsv"), data.table = FALSE)
mv <- mv[mv$tissue %in% c("Tumor", "Non-Tumor"), , drop = FALSE]
mv <- mv[mv$gsm %in% colnames(Xv), , drop = FALSE]
mv <- mv[!duplicated(mv$gsm), , drop = FALSE]
Xv <- Xv[, mv$gsm, drop = FALSE]
mv$tissue <- factor(mv$tissue, levels = c("Non-Tumor", "Tumor"))
mv$sex <- factor(mv$sex)
mv$age <- as.numeric(mv$age)
cat("GSE14520:", ncol(Xv), "amostras;", nrow(Xv), "genes;", sum(mv$tissue == "Tumor"), "tumores\n")

design <- model.matrix(~ 0 + tissue, data = mv)
colnames(design) <- c("NonTumor", "Tumor")
dup <- duplicateCorrelation(Xv, design, block = mv$pair)
fit <- lmFit(Xv, design, block = mv$pair, correlation = dup$consensus)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Tumor - NonTumor, levels = design)))
ttv <- topTable(fit2, coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
ttv$gene <- rownames(ttv)
ttv$regulation <- ifelse(ttv$adj.P.Val < FDR & ttv$logFC > LFC, "Up",
                  ifelse(ttv$adj.P.Val < FDR & ttv$logFC < -LFC, "Down", "NS"))
fwrite(ttv, file.path(dir_val, "GSE14520_tumor_vs_adjacente.csv"))
cat("GSE14520 DEGs tumor vs adjacente:", sum(ttv$regulation != "NS"), "\n")

# TCGA referência
ttc <- fread(file.path(dir_de, "modelo_tumor_vs_adjacente_base_completo.csv"), data.table = FALSE)
comuns <- intersect(ttc$gene, ttv$gene)
a <- deg_set(ttc); b <- deg_set(ttv)
ov <- intersect(a, b)
N <- length(comuns); n1 <- length(intersect(a, comuns)); n2 <- length(intersect(b, comuns))
p_hiper <- phyper(length(ov) - 1, n1, N - n1, n2, lower.tail = FALSE)
sp <- cor(ttc$logFC[match(comuns, ttc$gene)], ttv$logFC[match(comuns, ttv$gene)],
          method = "spearman", use = "complete.obs")
dir_ok <- mean(sign(ttc$logFC[match(ov, ttc$gene)]) == sign(ttv$logFC[match(ov, ttv$gene)]))
val_sig <- mean(ttv$regulation[match(a, ttv$gene)] != "NS", na.rm = TRUE)
val_same <- mean(sign(ttc$logFC[match(a, ttc$gene)]) == sign(ttv$logFC[match(a, ttv$gene)]), na.rm = TRUE)

val_tab <- data.frame(
  contraste = "Tumor vs adjacente",
  coorte_treinamento = "TCGA-LIHC (n = 418)",
  coorte_validacao = "GSE14520 (n = 445)",
  genes_comuns = N,
  degs_tcga = length(a), degs_validacao = length(b), degs_comuns = length(ov),
  jaccard = length(ov) / length(union(a, b)),
  hipergeometrico_p = p_hiper,
  spearman_logfc = sp,
  direcao_concordante = dir_ok,
  fracao_degs_tcga_significativos_na_validacao = val_sig,
  fracao_degs_tcga_mesma_direcao = val_same,
  stringsAsFactors = FALSE)
fwrite(val_tab, file.path(dir_val, "validacao_tumor_vs_adjacente.csv"))

rep <- merge(ttc[ttc$regulation != "NS", c("gene","logFC","adj.P.Val")],
             ttv[, c("gene","logFC","adj.P.Val","regulation")], by = "gene",
             suffixes = c("_tcga","_validacao"))
rep$mesma_direcao <- sign(rep$logFC_tcga) == sign(rep$logFC_validacao)
rep <- rep[order(rep$adj.P.Val_validacao), ]
fwrite(rep, file.path(dir_val, "GSE14520_genes_replicados.csv"))

# ---------------------------------------------------------------------------
# GSE14520: estágio (validação negativa esperada)
# ---------------------------------------------------------------------------
d_st <- mv[mv$tissue == "Tumor" & mv$stage_group %in% c("I-II", "III-IV"), , drop = FALSE]
d_st$stage_group <- factor(d_st$stage_group, levels = c("I-II", "III-IV"))
d_st <- d_st[complete.cases(d_st[, c("stage_group", "sex", "age")]), , drop = FALSE]
Xst <- Xv[, d_st$gsm, drop = FALSE]
covs_st <- if (length(unique(d_st$sex[!is.na(d_st$sex)])) > 1) c("sex", "age") else c("age")
design_st <- model.matrix(as.formula(paste("~ 0 + stage_group +", paste(covs_st, collapse = " + "))), data = d_st)
colnames(design_st) <- sub("stage_group", "group", colnames(design_st))
fit_st <- lmFit(Xst, design_st)
cst <- setNames(rep(0, ncol(design_st)), colnames(design_st))
cst["groupIII-IV"] <- 1; cst["groupI-II"] <- -1
tts <- topTable(eBayes(contrasts.fit(fit_st, cst), trend = TRUE, robust = TRUE),
                coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
tts$gene <- rownames(tts)
tts$regulation <- ifelse(tts$adj.P.Val < FDR & tts$logFC > LFC, "Up",
                  ifelse(tts$adj.P.Val < FDR & tts$logFC < -LFC, "Down", "NS"))
fwrite(tts, file.path(dir_val, "GSE14520_estagio.csv"))
cat("GSE14520 DEGs estágio:", sum(tts$regulation != "NS"), "de", nrow(d_st), "tumores\n")

ttc_st <- fread(file.path(dir_de, "modelo_estagio_base_completo.csv"), data.table = FALSE)
sig_st <- deg_set(ttc_st)
sig_up <- ttc_st$gene[ttc_st$regulation == "Up"]; sig_dn <- ttc_st$gene[ttc_st$regulation == "Down"]
up_v <- intersect(sig_up, rownames(Xst)); dn_v <- intersect(sig_dn, rownames(Xst))
Z <- t(scale(t(Xst)))
score_st <- colMeans(Z[up_v, , drop = FALSE], na.rm = TRUE) - colMeans(Z[dn_v, , drop = FALSE], na.rm = TRUE)
auc_st <- auc(score_st, d_st$stage_group == "III-IV")

com_st <- intersect(ttc_st$gene, tts$gene)
a_st <- intersect(sig_st, com_st); b_st <- intersect(deg_set(tts), com_st)
ov_st <- intersect(a_st, b_st)
p_st <- phyper(length(ov_st) - 1, length(a_st), length(com_st) - length(a_st),
               length(b_st), lower.tail = FALSE)
val_st <- data.frame(
  contraste = "Estágio III-IV vs I-II",
  coorte_treinamento = "TCGA-LIHC (n = 344)",
  coorte_validacao = sprintf("GSE14520 (n = %d)", nrow(d_st)),
  genes_comuns = length(com_st), degs_tcga = length(a_st), degs_validacao = length(b_st),
  degs_comuns = length(ov_st), jaccard = length(ov_st)/length(union(a_st, b_st)),
  hipergeometrico_p = p_st,
  spearman_logfc = cor(ttc_st$logFC[match(com_st, ttc_st$gene)],
                       tts$logFC[match(com_st, tts$gene)], method = "spearman", use = "complete.obs"),
  auc_assinatura_tcga = auc_st,
  genes_assinatura_usados = length(up_v) + length(dn_v),
  stringsAsFactors = FALSE)
fwrite(val_st, file.path(dir_val, "validacao_estagio.csv"))

# ---------------------------------------------------------------------------
# TCGA: validação cruzada do grau (5-fold, 20 repetições)
# ---------------------------------------------------------------------------
meta <- preparar_metadados(file.path(ROOT, "dados", "processed", "metadata_analise.tsv"))$meta
X <- carregar_expressao(file.path(ROOT, "dados", "processed", "expressao_log2tpm_filtrada.tsv.gz"))
tum <- meta[meta$group == "TCGA_Tumor" & !is.na(meta$grade_group), , drop = FALSE]
tum$grade_group <- factor(tum$grade_group, levels = c("G1-G2", "G3-G4"))
tum$age <- as.numeric(tum$age)

rodar_cv <- function(rep_i) {
  set.seed(2026 + rep_i)
  folds <- sample(rep(1:5, length.out = nrow(tum)))
  res <- lapply(1:5, function(k) {
    tr <- tum[folds != k, , drop = FALSE]; te <- tum[folds == k, , drop = FALSE]
    tt <- rodar_limma(X, tr, "grade_group", c("G1-G2", "G3-G4"), c("sex", "age"), fdr_cut = FDR, lfc_cut = LFC)
    if (is.null(tt)) return(data.frame(rep = rep_i, fold = k, auc = NA, n_deg = NA))
    up <- tt$gene[tt$regulation == "Up"]; dn <- tt$gene[tt$regulation == "Down"]
    if (length(up) + length(dn) == 0) return(data.frame(rep = rep_i, fold = k, auc = NA, n_deg = 0))
    mu <- rowMeans(X[, tr$sample, drop = FALSE]); sdv <- apply(X[, tr$sample, drop = FALSE], 1, sd)
    sdv[sdv == 0] <- 1
    Zt <- (X[, te$sample, drop = FALSE] - mu) / sdv
    score <- colMeans(Zt[up, , drop = FALSE], na.rm = TRUE) - colMeans(Zt[dn, , drop = FALSE], na.rm = TRUE)
    data.frame(rep = rep_i, fold = k, auc = auc(score, te$grade_group == "G3-G4"),
               n_deg = length(up) + length(dn))
  })
  do.call(rbind, res)
}
cv <- do.call(rbind, lapply(1:20, rodar_cv))
fwrite(cv, file.path(dir_val, "validacao_cruzada_grau.csv"))
cv_sum <- data.frame(
  repeticoes = 20, folds = 5,
  auc_media = mean(cv$auc, na.rm = TRUE),
  auc_ic95_inf = quantile(cv$auc, 0.025, na.rm = TRUE),
  auc_ic95_sup = quantile(cv$auc, 0.975, na.rm = TRUE),
  auc_min = min(cv$auc, na.rm = TRUE), auc_max = max(cv$auc, na.rm = TRUE),
  folds_sem_deg = sum(cv$n_deg == 0, na.rm = TRUE), folds_total = nrow(cv))
fwrite(cv_sum, file.path(dir_val, "validacao_cruzada_grau_resumo.csv"))

centros <- names(which(table(tum$tss_collapsed) >= 20))
loco <- lapply(centros, function(cc) {
  tr <- tum[tum$tss_collapsed != cc, , drop = FALSE]; te <- tum[tum$tss_collapsed == cc, , drop = FALSE]
  tt <- rodar_limma(X, tr, "grade_group", c("G1-G2", "G3-G4"), c("sex", "age"), fdr_cut = FDR, lfc_cut = LFC)
  if (is.null(tt)) return(data.frame(centro = cc, n_treino = nrow(tr), n_teste = nrow(te), auc = NA, n_deg = NA))
  up <- tt$gene[tt$regulation == "Up"]; dn <- tt$gene[tt$regulation == "Down"]
  if (length(up) + length(dn) == 0) return(data.frame(centro = cc, n_treino = nrow(tr), n_teste = nrow(te), auc = NA, n_deg = 0))
  mu <- rowMeans(X[, tr$sample, drop = FALSE]); sdv <- apply(X[, tr$sample, drop = FALSE], 1, sd); sdv[sdv == 0] <- 1
  Zt <- (X[, te$sample, drop = FALSE] - mu) / sdv
  score <- colMeans(Zt[up, , drop = FALSE], na.rm = TRUE) - colMeans(Zt[dn, , drop = FALSE], na.rm = TRUE)
  data.frame(centro = cc, n_treino = nrow(tr), n_teste = nrow(te),
             auc = auc(score, te$grade_group == "G3-G4"), n_deg = length(up) + length(dn))
})
loco <- do.call(rbind, loco)
fwrite(loco, file.path(dir_val, "validacao_leave_one_center_out_grau.csv"))
cat("\nResumo CV grau:\n"); print(cv_sum)
cat("\nLOCO grau:\n"); print(loco)
