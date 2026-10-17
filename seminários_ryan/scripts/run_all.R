# =============================================================================
# run_all.R
# -----------------------------------------------------------------------------
# Executa todo o pipeline do seminário, em ordem.
#
# Uso (a partir de seminários_ryan/):
#   Rscript scripts/run_all.R
#
# Para usar o transcriptoma completo, defina antes:
#   Sys.setenv(TCGA_GTEX_TPM = "C:/caminho/TcgaTargetGtex_rsem_gene_tpm.gz")
# ou exporte a variável de ambiente no terminal.
#
# Se o transcriptoma completo não estiver disponível, a etapa 00 usa o painel
# LA (213 genes) e as etapas seguintes rodam com esse painel. A etapa 04
# (demonstração do ComBat) só roda se a matriz completa existir.
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
RSCRIPT <- file.path(R.home("bin"), "Rscript")

rodar <- function(script) {
  caminho <- file.path(SCRIPT_DIR, script)
  if (!file.exists(caminho)) stop("Script não encontrado: ", caminho)
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(">> ", script, "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")
  status <- system2(RSCRIPT, c("--vanilla", shQuote(caminho)))
  if (!identical(status, 0L)) stop("Falha em ", script, " (status ", status, ")")
  invisible(TRUE)
}

# Caminho do TPM completo: variável de ambiente ou candidato no TEMP
if (!nzchar(Sys.getenv("TCGA_GTEX_TPM"))) {
  candidatos <- c(
    file.path(Sys.getenv("TEMP"), "hepato_full", "TcgaTargetGtex_rsem_gene_tpm.gz"),
    file.path(ROOT, "dados", "raw", "TcgaTargetGtex_rsem_gene_tpm.gz")
  )
  candidatos <- candidatos[file.exists(candidatos)]
  if (length(candidatos) > 0) Sys.setenv(TCGA_GTEX_TPM = candidatos[1])
}

cat("Raiz:", ROOT, "\n")
cat("TPM completo:", ifelse(nzchar(Sys.getenv("TCGA_GTEX_TPM")),
                             Sys.getenv("TCGA_GTEX_TPM"), "(não encontrado)"), "\n")

rodar("00_preparar_dados.R")
rodar("01_analise_exploratoria.R")
rodar("02_expressao_diferencial.R")
rodar("03_volcano_plots.R")

mat_completa <- file.path(ROOT, "dados", "processed", "expressao_log2tpm_filtrada.tsv.gz")
if (file.exists(mat_completa)) {
  rodar("04_ComBat_demonstracao.R")
} else {
  cat("\n[SKIP] 04_ComBat_demonstracao.R: matriz completa não encontrada.\n")
}

rodar("05_fenotipos_avancados.R")
rodar("10_auditoria_dados_e_volcano.R")
rodar("11_benchmarks.R")
rodar("12_compliance.R")

cat("\nPipeline concluído. Resultados em:", file.path(ROOT, "resultados"), "\n")
cat("Auditorias em:", file.path(ROOT, "auditorias"), "\n")
