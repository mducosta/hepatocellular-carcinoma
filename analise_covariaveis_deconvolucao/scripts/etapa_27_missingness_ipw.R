#!/usr/bin/env Rscript
# =============================================================================
# etapa_27_missingness_ipw.R
# Sensibilidade à ausência de dados por ponderação de probabilidade inversa.
#
# Para contrastes cuja variável de grupo tem missingness alta (sorologia,
# Child-Pugh, fibrose e invasão vascular), o modelo de casos completos pode
# representar uma subpopulação selecionada. Aqui:
#   1. modela-se a probabilidade de observação do grupo com variáveis auxiliares
#      completas (sexo, idade, estágio, grau, histologia e TSS);
#   2. calcula-se peso estabilizado w = P(observado)/P(observado | aux);
#   3. repete-se o limma com pesos amostrais;
#   4. compara-se o número de DEGs e a ESS (effective sample size).
# A ponderação não cria informação: ela reequilibra os casos observados para a
# população auxiliar. Se o contraste continuar nulo, o resultado é mais crível.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(limma); library(ggplot2)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))

dir_proc <- file.path(ROOT, "dados", "processed")
dir_out  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "sensibilidades")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
FDR_CUT <- 0.05; LFC_CUT <- 1

prep <- preparar_metadados(file.path(dir_proc, "metadata_analise.tsv"))
meta <- prep$meta
X <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))
tum <- meta[meta$group == "TCGA_Tumor", , drop = FALSE]
tum$age <- as.numeric(tum$age)
tum$histology_group <- factor(ifelse(is.na(tum$histology_group), "HCC", as.character(tum$histology_group)),
                              levels = c("HCC", "Other"))
tum$stage_group <- factor(ifelse(is.na(tum$stage_group), NA, as.character(tum$stage_group)),
                          levels = c("I-II", "III-IV"))
tum$grade_group <- factor(ifelse(is.na(tum$grade_group), NA, as.character(tum$grade_group)),
                          levels = c("G1-G2", "G3-G4"))
aux_vars <- c("sex", "age", "stage_group", "grade_group", "histology_group", "tss_collapsed")
aux_complete <- complete.cases(tum[, aux_vars])
cat("Tumores com auxiliares completos:", sum(aux_complete), "de", nrow(tum), "\n")

# Contraste: coluna de grupo, níveis, rótulo
contrastes <- list(
  list(nome = "estagio", col = "stage_group", níveis = c("I-II", "III-IV")),
  list(nome = "grau", col = "grade_group", níveis = c("G1-G2", "G3-G4")),
  list(nome = "sorologia_hbv_vs_hcv", col = "viral_known", níveis = c("HBV", "HCV")),
  list(nome = "child_pugh", col = "child_pugh_group", níveis = c("A", "B-C")),
  list(nome = "fibrose", col = "fibrosis_adv2", níveis = c("Ishak < 5", "Ishak >= 5")),
  list(nome = "invasao_vascular", col = "vascular_invasion2", níveis = c("Micro", "Macro"))
)

fit_weighted <- function(Xsub, d, grupo, níveis, covs, w) {
  dd <- d[d[[grupo]] %in% níveis, , drop = FALSE]
  dd[[grupo]] <- factor(as.character(dd[[grupo]]), levels = níveis)
  keep <- complete.cases(dd[, c(grupo, covs), drop = FALSE])
  ww <- w[match(dd$sample, d$sample)]
  dd <- dd[keep, , drop = FALSE]; ww <- ww[keep]
  if (nrow(dd) < 8) return(NULL)
  for (v in names(dd)) if (is.factor(dd[[v]])) dd[[v]] <- droplevels(dd[[v]])
  fml <- as.formula(paste("~ 0 +", grupo, "+", paste(covs, collapse = " + ")))
  design <- model.matrix(fml, data = dd)
  colnames(design) <- sub(grupo, "group", colnames(design), fixed = TRUE)
  if (qr(design)$rank < ncol(design)) return(NULL)
  W <- matrix(ww, nrow = nrow(Xsub), ncol = length(ww), byrow = TRUE)
  fit <- lmFit(Xsub[, dd$sample, drop = FALSE], design, weights = W)
  cvec <- setNames(rep(0, ncol(design)), colnames(design))
  cvec[paste0("group", níveis[1])] <- -1
  cvec[paste0("group", níveis[2])] <- 1
  fit2 <- eBayes(contrasts.fit(fit, cvec), trend = TRUE, robust = TRUE)
  tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt$regulation <- ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC > LFC_CUT, "Up",
                   ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC < -LFC_CUT, "Down", "NS"))
  tt
}

res <- list()
for (cs in contrastes) {
  d <- tum
  d$.target <- d[[cs$col]]
  d$.observed <- !is.na(d$.target)
  # Modelo de observação com auxiliares completos
  d_aux <- d[aux_complete, , drop = FALSE]
  d_aux$.observed <- as.integer(!is.na(d_aux[[cs$col]]))
  if (length(unique(d_aux$.observed)) < 2) {
    cat("Contraste", cs$nome, "sem variação de observação; ignorado.\n"); next
  }
  fml_obs <- as.formula(paste(".observed ~", paste(aux_vars, collapse = " + ")))
  mod_obs <- glm(fml_obs, data = d_aux, family = binomial())
  p_obs <- predict(mod_obs, newdata = d, type = "response")
  p_obs <- pmin(pmax(p_obs, 0.02), 0.98)
  p_bar <- mean(d$.observed, na.rm = TRUE)
  w <- ifelse(d$.observed, p_bar / p_obs, NA_real_)
  w <- pmin(w, quantile(w, 0.99, na.rm = TRUE))  # evita pesos extremos
  ess <- sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)

  base <- rodar_limma(X, d, cs$col, cs$níveis, c("sex", "age"),
                      fdr_cut = FDR_CUT, lfc_cut = LFC_CUT)
  wfit <- fit_weighted(X, d, cs$col, cs$níveis, c("sex", "age"), w)
  n_complete <- sum(!is.na(d[[cs$col]]) & complete.cases(d[, c(cs$col, "sex", "age")]))
  n_missing <- sum(is.na(d[[cs$col]]))
  deg_b <- if (is.null(base)) NA else sum(base$regulation != "NS")
  deg_w <- if (is.null(wfit)) NA else sum(wfit$regulation != "NS")
  jac <- if (!is.null(base) && !is.null(wfit)) {
    a <- base$gene[base$regulation != "NS"]; b <- wfit$gene[wfit$regulation != "NS"]
    if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA
  } else NA
  res[[cs$nome]] <- data.frame(
    contraste = cs$nome, n_complete = n_complete, n_missing = n_missing,
    pct_missing = round(100 * n_missing / nrow(d), 1),
    n_aux_completo = nrow(d_aux),
    auc_modelo_observacao = {pr <- fitted(mod_obs); r <- rank(pr); n1 <- sum(d_aux$.observed == 1); n0 <- sum(d_aux$.observed == 0); as.numeric((sum(r[d_aux$.observed == 1]) - n1*(n1+1)/2)/(n1*n0))},
    peso_min = min(w, na.rm = TRUE), peso_max = max(w, na.rm = TRUE),
    ess = round(ess, 1), degs_sem_peso = deg_b, degs_com_peso = deg_w,
    jaccard = jac, stringsAsFactors = FALSE)
  if (!is.null(wfit)) fwrite(wfit, file.path(dir_out, paste0("ipw_", cs$nome, "_completo.csv")))
  cat("Contraste", cs$nome, "ok: sem peso", deg_b, "com peso", deg_w, "\n")
}
resumo <- do.call(rbind, res)
fwrite(resumo, file.path(dir_out, "sensibilidade_missingness_ipw.csv"))
cat("\nResumo IPW:\n"); print(resumo, row.names = FALSE)
