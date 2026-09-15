# =============================================================================
# Nanotec 2026 — GINAB/UENF
# Etapa 03 — Sobrevivência (Kaplan-Meier + Cox) para MMP9 e MMPs relevantes
# -----------------------------------------------------------------------------
# Fonte de sobrevivência: TCGA_survival_data (colunas OS / OS.time / DSS / PFI).
# Amostras: apenas o grupo Tumor (TCGA-LIHC Primary Tumor, n=369).
# Estratificação: mediana da expressão (log2) de cada gene no grupo tumoral.
# Saídas (outputs/survival/):
#   KM_OS_<gene>.png, survival_logrank.csv, cox_univariado.csv
#
# IMPORTANTE: trata-se de ASSOCIAÇÃO, não causalidade.
#
# Executar a partir da pasta nanotec2026/:
#   Rscript scripts/03_survival.R
# =============================================================================

cat("============================================================\n")
cat("NANOTEC 2026 - Etapa 03: Sobrevivência (KM + Cox)\n")
cat("============================================================\n")

base     <- "."
dir_proc <- file.path(base, "dados", "processed")
dir_raw  <- file.path(base, "dados", "raw")
dir_surv <- file.path(base, "outputs", "survival")
dir.create(dir_surv, showWarnings = FALSE, recursive = TRUE)

for (p in c("survival", "survminer")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))

# -----------------------------------------------------------------------------
# 1. Leitura
# -----------------------------------------------------------------------------
mat <- read.delim(file.path(dir_proc, "liver_panel_41genes.tsv"),
                  header = TRUE, sep = "\t", check.names = FALSE, quote = "")
rownames(mat) <- mat[[1]]; mat <- mat[, -1, drop = FALSE]
annot <- read.delim(file.path(dir_proc, "samples_annotation.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE, quote = "")
m <- match(colnames(mat), annot[["sample"]])
grupo <- annot[["grupo"]][m]

# survival data (TCGA apenas)
f_surv <- list.files(dir_raw, pattern = "TCGA_survival_data", full.names = TRUE)[1]
surv <- read.delim(f_surv, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
cat(sprintf("[SURV] %d amostras no arquivo de sobrevivência.\n", nrow(surv)))

# amostras tumorais com dado de sobrevivência
tumor_samples <- colnames(mat)[grupo == "Tumor"]
idx <- match(tumor_samples, surv[["sample"]])
tem_surv <- !is.na(idx) & !is.na(surv[["OS.time"]][idx]) & !is.na(surv[["OS"]][idx])
cat(sprintf("[SURV] %d tumores LIHC com OS/OS.time disponíveis (de %d).\n",
            sum(tem_surv), length(tumor_samples)))

expr_tumor <- mat[, tumor_samples[tem_surv], drop = FALSE]
os_time <- surv[["OS.time"]][idx[tem_surv]]
os_evt  <- surv[["OS"]][idx[tem_surv]]

# -----------------------------------------------------------------------------
# 2. Kaplan-Meier + log-rank + Cox (univariado) por gene
# -----------------------------------------------------------------------------
genes_alvo <- c("MMP9", "MMP1", "MMP11", "MMP12", "MMP14")  # MMPs up significativos
genes_alvo <- genes_alvo[genes_alvo %in% rownames(expr_tumor)]

res_logrank <- data.frame(gene = character(), p_logrank = numeric(),
                          n_high = integer(), n_low = integer(),
                          stringsAsFactors = FALSE)
res_cox <- data.frame(gene = character(), HR = numeric(), HR_low = numeric(),
                      HR_high = numeric(), p = numeric(), stringsAsFactors = FALSE)

for (g in genes_alvo) {
  x <- as.numeric(expr_tumor[g, ])
  cutoff <- median(x, na.rm = TRUE)
  grp_hl <- factor(ifelse(x >= cutoff, "High", "Low"), levels = c("Low", "High"))

  # Kaplan-Meier
  sfit <- survfit(Surv(os_time, os_evt) ~ grp_hl)
  lr <- survdiff(Surv(os_time, os_evt) ~ grp_hl)
  p_lr <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

  png(file.path(dir_surv, paste0("KM_OS_", g, ".png")), width = 900, height = 700, res = 130)
  plt <- ggsurvplot(sfit, data = data.frame(os_time, os_evt, grp_hl),
                    pval = TRUE, risk.table = TRUE,
                    title = paste0("Sobrevida global — ", g, " (LIHC, mediana)"),
                    xlab = "Tempo (dias)", ylab = "Probabilidade de sobrevida",
                    legend.title = g, palette = c("#4575b4", "#d73027"))
  print(plt)
  dev.off()

  # Cox univariado (expressão contínua, por unidade de log2)
  cox <- coxph(Surv(os_time, os_evt) ~ x)
  ci <- exp(confint(cox))
  res_cox <- rbind(res_cox, data.frame(gene = g, HR = exp(coef(cox)),
                                       HR_low = ci[1], HR_high = ci[2],
                                       p = summary(cox)$coefficients[1, "Pr(>|z|)"],
                                       stringsAsFactors = FALSE))
  res_logrank <- rbind(res_logrank, data.frame(gene = g, p_logrank = p_lr,
                                               n_high = sum(grp_hl == "High"),
                                               n_low = sum(grp_hl == "Low"),
                                               stringsAsFactors = FALSE))
  cat(sprintf("[KM/Cox] %s: log-rank p=%.4g, Cox HR=%.3f (%.3f-%.3f) p=%.4g\n",
              g, p_lr, exp(coef(cox)), ci[1], ci[2],
              summary(cox)$coefficients[1, "Pr(>|z|)"]))
}

write.csv(res_logrank, file.path(dir_surv, "survival_logrank.csv"), row.names = FALSE)
write.csv(res_cox, file.path(dir_surv, "cox_univariado.csv"), row.names = FALSE)

cat("\n[FIM] Etapa 03 concluída. Resultados em outputs/survival/.\n")
