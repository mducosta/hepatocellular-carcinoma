# =============================================================================
# etapa_00_baixar_dados_brutos.R
# -----------------------------------------------------------------------------
# Download dos dados brutos usados no seminário.
#
# Por padrão, baixa apenas os arquivos pequenos (clínico/fenótipo). O
# transcriptoma completo do Toil (~1,3 GB) só é baixado se a variável de
# ambiente BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE.
#
# Uso:
#   Rscript scripts/etapa_00_baixar_dados_brutos.R
#   BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE Rscript scripts/etapa_00_baixar_dados_brutos.R
#
# Depois do download, rode:
#   TCGA_GTEX_TPM="<caminho>/TcgaTargetGtex_rsem_gene_tpm.gz" \
#     Rscript scripts/etapa_01_preparar_matrizes_e_metadados.R
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_raw <- file.path(ROOT, "dados", "raw")
dir.create(dir_raw, recursive = TRUE, showWarnings = FALSE)

baixar <- function(url, destino, obrigatorio = TRUE) {
  if (file.exists(destino) && file.info(destino)$size > 0) {
    cat("[SKIP] Já existe:", destino, "\n")
    return(invisible(TRUE))
  }
  cat("[DOWNLOAD]", url, "\n")
  ok <- tryCatch({
    utils::download.file(url, destino, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) {
    warning("Falha ao baixar ", url, ": ", e$message)
    FALSE
  })
  if (!ok && obrigatorio) stop("Download obrigatório falhou: ", url)
  invisible(ok)
}

# -----------------------------------------------------------------------------
# 1. Arquivos pequenos (sempre)
# -----------------------------------------------------------------------------
baixar(
  "https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/LIHC_clinicalMatrix",
  file.path(dir_raw, "LIHC_clinicalMatrix")
)

baixar(
  "https://gdc.xenahubs.net/download/TCGA-LIHC.clinical.tsv.gz",
  file.path(dir_raw, "TCGA-LIHC.clinical.tsv.gz")
)

baixar(
  "https://toil.xenahubs.net/download/GTEX_phenotype.gz",
  file.path(dir_raw, "GTEX_phenotype.gz")
)

baixar(
  "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz",
  file.path(dir_raw, "TcgaTargetGTEX_phenotype.txt.gz")
)

baixar(
  "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
  file.path(dir_raw, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
)

# -----------------------------------------------------------------------------
# 2. Transcriptoma completo (opcional; ~1,3 GB)
# -----------------------------------------------------------------------------
if (toupper(Sys.getenv("BAIXAR_TRANSCRIPTOMA_COMPLETO", "FALSE")) %in% c("TRUE", "1", "SIM", "YES")) {
  baixar(
    "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    file.path(dir_raw, "TcgaTargetGtex_rsem_gene_tpm.gz")
  )
} else {
  cat("\n[INFO] Transcriptoma completo não baixado (1,3 GB).\n")
  cat("       Defina BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE para baixar.\n")
  cat("       Alternativa: use a matriz processada já incluída em\n")
  cat("       dados/processed/expressao_log2tpm_filtrada.tsv.gz.\n")
}

cat("\nDownloads concluídos em:", dir_raw, "\n")
