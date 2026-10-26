# =============================================================================
# executar_todas_as_etapas.R
# -----------------------------------------------------------------------------
# Executa o pipeline completo do seminário, em ordem.
#
# Uso (a partir de seminários_ryan/):
#   Rscript scripts/executar_todas_as_etapas.R
#
# As etapas 13–19 correspondem às análises avançadas (covariáveis, ComBat,
# deconvolução, integração, auditoria, benchmark e compliance). A geração do
# artigo final não faz parte deste pipeline público.
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

PYTHON <- Sys.which(c("python", "python3"))[1]
rodar_py <- function(script) {
  caminho <- file.path(SCRIPT_DIR, script)
  if (!file.exists(caminho)) stop("Script não encontrado: ", caminho)
  cat("
", strrep("=", 78), "
", sep = "")
  cat(">> ", script, "
", sep = "")
  cat(strrep("=", 78), "
", sep = "")
  status <- system2(PYTHON, c(shQuote(caminho)))
  if (!identical(status, 0L)) stop("Falha em ", script, " (status ", status, ")")
  invisible(TRUE)
}

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

# Pipeline exploratório e clássico
rodar("etapa_01_preparar_matrizes_e_metadados.R")
rodar("etapa_02_analise_exploratoria_fenotipos.R")
rodar("etapa_03_expressao_diferencial_limma.R")
rodar("etapa_04_volcano_plots.R")

mat_completa <- file.path(ROOT, "dados", "processed", "expressao_log2tpm_filtrada.tsv.gz")
if (file.exists(mat_completa)) {
  rodar("etapa_05_correcao_lote_ComBat.R")
} else {
  cat("\n[SKIP] etapa_05_correcao_lote_ComBat.R: matriz completa não encontrada.\n")
}

rodar("etapa_06_fenotipos_clinicos_avancados.R")
rodar("etapa_07_grupos_alternativos.R")
rodar("etapa_08_analises_grupos_separados.R")
rodar("etapa_10_auditoria_dados_volcano.R")
rodar("etapa_11_benchmarks.R")
rodar("etapa_12_compliance.R")

# Análises avançadas
rodar("etapa_13_expressao_diferencial_covariaveis.R")
rodar("etapa_14_ComBat_testes_multiplos.R")
rodar("etapa_15_deconvolucao_celular.R")
rodar("etapa_16_integracao_deconvolucao_expressao.R")
rodar("etapa_16b_via_LA_covariaveis.R")
rodar_py("etapa_17_auditoria_analises_avancadas.py")
rodar("etapa_18_benchmarks_analises_avancadas.R")
rodar_py("etapa_19_compliance_analises_avancadas.py")

cat("\nPipeline concluído. Resultados em:", file.path(ROOT, "resultados"), "\n")
cat("Auditorias em:", file.path(ROOT, "auditorias"), "\n")
cat("Artigo final em:", file.path(ROOT, "artigo_final"), "\n")
