# =============================================================================
# Nanotec 2026 — GINAB/UENF
# ORQUESTRADOR — roda as etapas na ordem correta
# -----------------------------------------------------------------------------
# Uso (a partir da pasta nanotec2026/):
#   Rscript scripts/run_all.R
#
# Etapas:
#   00 - Download (UCSCXenaTools) + extração + mapeamento Ensembl->símbolo
#   01 - Auditoria (para se houver dados insuficientes)
#   02 - Expressão diferencial + volcano + heatmap + correlação MMP9
#   03 - Sobrevivência (KM + Cox)
#   04 - GSVA/ssGSEA de vias ECM
# =============================================================================

base <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))[1]], winslash = "/"))
base <- dirname(base)  # scripts/ -> nanotec2026/
setwd(base)
cat("Diretório de trabalho:", getwd(), "\n\n")

scripts <- c(
  "scripts/00_download_data.R",
  "scripts/01_audit_data.R",
  "scripts/02_de_analysis.R",
  "scripts/03_survival.R",
  "scripts/04_enrichment.R"
)

for (s in scripts) {
  cat("\n############################################################\n")
  cat("## EXECUTANDO:", s, "\n")
  cat("############################################################\n\n")
  status <- system2(file.path(R.home("bin"), "Rscript"), s, wait = TRUE)
  if (status != 0) {
    stop("Falha na etapa ", s, " (status ", status, "). Pipeline interrompido.")
  }
}

cat("\n============================================================\n")
cat("PIPELINE NANOTEC 2026 CONCLUÍDO COM SUCESSO.\n")
cat("============================================================\n")
