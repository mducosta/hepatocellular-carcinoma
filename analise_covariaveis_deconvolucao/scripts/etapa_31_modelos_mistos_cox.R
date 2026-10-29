#!/usr/bin/env Rscript
# =============================================================================
# etapa_31_modelos_mistos_cox.R
# (i) Modelos mistos com TSS como efeito aleatório para estágio e grau;
# (ii) recálculo dos modelos de Cox com casos completos e relato de cox.zph.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(lmerTest); library(survival); library(parallel)
})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))
dir_out <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "sensibilidades")
dir_tab <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
FDR <- 0.05; LFC <- 1

prep <- preparar_metadados(file.path(ROOT, "dados", "processed", "metadata_analise.tsv"))
meta <- prep$meta
X <- carregar_expressao(file.path(ROOT, "dados", "processed", "expressao_log2tpm_filtrada.tsv.gz"))
tum <- meta[meta$group == "TCGA_Tumor", , drop = FALSE]
tum$age <- as.numeric(tum$age)

# ---------------------------------------------------------------------------
# Modelos mistos: duplicateCorrelation (todos os genes) + lmerTest (top genes)
# ---------------------------------------------------------------------------
rodar_misto <- function(grupo, niveis, nome) {
  d <- tum[tum[[grupo]] %in% niveis & !is.na(tum[[grupo]]), , drop = FALSE]
  d[[grupo]] <- factor(as.character(d[[grupo]]), levels = niveis)
  d <- d[complete.cases(d[, c(grupo, "sex", "age", "tss_collapsed")]), , drop = FALSE]
  design <- model.matrix(as.formula(paste("~ 0 +", grupo, "+ sex + age")), data = d)
  colnames(design) <- sub(grupo, "group", colnames(design), fixed = TRUE)
  # duplicateCorrelation: intercepto aleatório por TSS (estrutura compoundsymmetry)
  dup <- duplicateCorrelation(X[, d$sample], design, block = d$tss_collapsed)
  fit <- lmFit(X[, d$sample], design, block = d$tss_collapsed, correlation = dup$consensus)
  cvec <- setNames(rep(0, ncol(design)), colnames(design))
  cvec[paste0("group", niveis[1])] <- -1
  cvec[paste0("group", niveis[2])] <- 1
  tt <- topTable(eBayes(contrasts.fit(fit, cvec), trend = TRUE, robust = TRUE),
                 coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt$regulation <- ifelse(tt$adj.P.Val < FDR & tt$logFC > LFC, "Up",
                   ifelse(tt$adj.P.Val < FDR & tt$logFC < -LFC, "Down", "NS"))
  fwrite(tt, file.path(dir_out, paste0("misto_duplicateCorrelation_", nome, ".csv")))

  # lmerTest nos 300 genes mais significativos do modelo fixo
  top <- head(tt$gene[order(tt$P.Value)], 300)
  dd <- d
  res_lmer <- lapply(top, function(g) {
    dd$.y <- as.numeric(X[g, dd$sample])
    m <- tryCatch(lmerTest::lmer(as.formula(paste(".y ~", grupo, "+ sex + age + (1 | tss_collapsed)")), data = dd),
                  error = function(e) NULL)
    if (is.null(m)) return(data.frame(gene = g, p_lmer = NA, logFC_lmer = NA))
    an <- tryCatch(anova(m, type = 3), error = function(e) NULL)
    if (is.null(an) || !("Pr(>F)" %in% names(an))) return(data.frame(gene = g, p_lmer = NA, logFC_lmer = NA))
    data.frame(gene = g, p_lmer = an[["Pr(>F)"]][1],
               logFC_lmer = fixef(m)[[paste0(grupo, niveis[2])]])
  })
  lmer_df <- do.call(rbind, res_lmer)
  lmer_df$p_lmer_bh <- p.adjust(lmer_df$p_lmer, method = "BH")
  lmer_df$p_limma <- tt$P.Value[match(lmer_df$gene, tt$gene)]
  lmer_df$logFC_limma <- tt$logFC[match(lmer_df$gene, tt$gene)]
  fwrite(lmer_df, file.path(dir_out, paste0("misto_lmer_", nome, ".csv")))
  data.frame(contraste = nome, n = nrow(d), genes = nrow(tt),
             consensus_correlation = dup$consensus,
             degs_duplicateCorrelation = sum(tt$regulation != "NS"),
             up = sum(tt$regulation == "Up"), down = sum(tt$regulation == "Down"),
             lmer_genes = nrow(lmer_df),
             lmer_degs_bh = sum(lmer_df$p_lmer_bh < FDR, na.rm = TRUE),
             spearman_p_lmer_vs_limma = cor(lmer_df$p_lmer, lmer_df$p_limma,
                                            method = "spearman", use = "complete.obs"),
             spearman_logfc_lmer_vs_limma = cor(lmer_df$logFC_lmer, lmer_df$logFC_limma,
                                                method = "spearman", use = "complete.obs"),
             stringsAsFactors = FALSE)
}
res_misto <- rbind(
  rodar_misto("stage_group", c("I-II", "III-IV"), "estagio"),
  rodar_misto("grade_group", c("G1-G2", "G3-G4"), "grau"))
fwrite(res_misto, file.path(dir_out, "modelos_mistos_resumo.csv"))
print(res_misto)

# ---------------------------------------------------------------------------
# Cox: casos completos, tempo > 0, cox.zph
# ---------------------------------------------------------------------------
surv_df <- tum[!is.na(tum$OS_time) & tum$OS_time > 0 & !is.na(tum$OS), , drop = FALSE]
surv_df$OS <- as.integer(surv_df$OS)
vars <- c("OS_time", "OS", "stage_group", "grade_group", "sex", "age")
dA <- surv_df[complete.cases(surv_df[, vars]), , drop = FALSE]
fitA <- coxph(Surv(OS_time, OS) ~ stage_group + grade_group + sex + age, data = dA)
tidy <- broom::tidy(fitA, exponentiate = TRUE, conf.int = TRUE)
tidy$n <- nrow(dA); tidy$eventos <- sum(dA$OS)
tidy$modelo <- "clinico_casos_completos"
fwrite(tidy, file.path(dir_tab, "cox_clinico_corrigido.csv"))

# modelo B com AFP e fibrose
dA$afp_high <- factor(ifelse(is.na(dA$afp), NA, ifelse(dA$afp >= 400, ">=400", "<400")),
                      levels = c("<400", ">=400"))
dA$fibrose <- factor(ifelse(is.na(dA$fibrosis_ishak), NA, ifelse(dA$fibrosis_ishak >= 5, ">=5", "<5")),
                     levels = c("<5", ">=5"))
dA$viral_known <- factor(ifelse(is.na(dA$viral_status) | dA$viral_status == "Unknown", NA,
                                as.character(dA$viral_status)))
varsB <- c(vars, "afp_high", "fibrose")
dB <- dA[complete.cases(dA[, varsB]), , drop = FALSE]
fitB <- coxph(Surv(OS_time, OS) ~ stage_group + grade_group + sex + age + afp_high + fibrose, data = dB)
tidyB <- broom::tidy(fitB, exponentiate = TRUE, conf.int = TRUE)
tidyB$n <- nrow(dB); tidyB$eventos <- sum(dB$OS); tidyB$modelo <- "clinico_afp_fibrose"
fwrite(tidyB, file.path(dir_tab, "cox_clinico_afp_fibrose.csv"))

# cox.zph
zphA <- as.data.frame(cox.zph(fitA)$table); zphA$modelo <- "clinico_casos_completos"; zphA$termo <- rownames(zphA)
zphB <- as.data.frame(cox.zph(fitB)$table); zphB$modelo <- "clinico_afp_fibrose"; zphB$termo <- rownames(zphB)
zph <- rbind(zphA, zphB)
fwrite(zph, file.path(dir_tab, "cox_zph.csv"))

# log-rank do estágio nos casos completos
sd <- survdiff(Surv(OS_time, OS) ~ stage_group, data = dA)
logrank <- data.frame(n = nrow(dA), eventos = sum(dA$OS), chisq = sd$chisq,
                      df = length(sd$n) - 1,
                      p = 1 - pchisq(sd$chisq, df = length(sd$n) - 1))
fwrite(logrank, file.path(dir_tab, "logrank_estagio.csv"))

# Cox da composição + cox.zph
uni <- fread(file.path(dir_tab, "deconvolucao", "deconvolucao_unificada_amostras.csv"),
             data.table = FALSE)
uni <- uni[uni$sample %in% surv_df$sample, , drop = FALSE]
uni$z_hep <- as.numeric(scale(uni[["xCell_Hepatocytes"]]))
cox_cel <- merge(surv_df[, c("sample", "OS", "OS_time", "age", "sex", "stage_group", "grade_group")],
                 uni[, c("sample", "z_hep")], by = "sample")
cox_cel <- cox_cel[complete.cases(cox_cel), , drop = FALSE]
fitC <- coxph(Surv(OS_time, OS) ~ z_hep + age + sex + stage_group + grade_group, data = cox_cel)
zc <- as.data.frame(cox.zph(fitC)$table); zc$modelo <- "cox_celulas_hepatocytes"; zc$termo <- rownames(zc)
fwrite(zc, file.path(dir_tab, "cox_zph_celulas.csv"))
fwrite(data.frame(modelo = "cox_celulas_hepatocytes", n = nrow(cox_cel), eventos = sum(cox_cel$OS),
                  termo = "z_hep", HR = exp(coef(fitC))[["z_hep"]],
                  IC95_inf = exp(confint(fitC))[["z_hep", 1]],
                  IC95_sup = exp(confint(fitC))[["z_hep", 2]],
                  p = summary(fitC)$coefficients["z_hep", "Pr(>|z|)"]),
       file.path(dir_tab, "cox_celulas_hepatocytes_corrigido.csv"))

cat("\nCox casos completos:\n"); print(tidy[, c("term", "estimate", "conf.low", "conf.high", "p.value", "n", "eventos")])
cat("\ncox.zph global:\n"); print(zph[zph$termo == "GLOBAL", c("modelo", "p")])
cat("\nlog-rank:\n"); print(logrank)
cat("\nCox células:\n"); print(zc)
