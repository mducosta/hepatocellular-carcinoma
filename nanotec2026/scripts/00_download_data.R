# =============================================================================
# Nanotec 2026 — I Workshop Interdisciplinar de Nanotecnologia e Inovação
# GINAB/UENF
#
# Etapa 00 — Download (UCSCXenaTools) + extração + mapeamento Ensembl->símbolo
# -----------------------------------------------------------------------------
# Tema: Remodelamento da MEC e expressão de MMP9/MMPs/TIMPs no HCC
# Fonte: UCSC Xena — hub Toil (https://toil.xenahubs.net)
#
# Datasets baixados:
#   1. TcgaTargetGtex_rsem_gene_tpm   -> expressão TCGA+GTEx (TPM; linhas = ENSG)
#   2. TcgaTargetGTEX_phenotype.txt   -> fenótipo (tecido / tipo de amostra)
#   3. TCGA_survival_data             -> sobrevivência (OS/OS.time, etc.)
#
# ETAPAS:
#   A. Baixa os 3 datasets com UCSCXenaTools (XenaGenerate->XenaQuery->XenaDownload).
#   B. Lê o fenótipo e seleciona amostras de fígado:
#        - TCGA-LIHC Primary Tumor   (Tumor)     n=369
#        - GTEx Liver Normal Tissue  (Normal)    n=110
#        - TCGA Solid Tissue Normal  (Adjacent)  n=50
#      (TCGA "Recurrent Tumor" n=2 é EXCLUÍDO: recidiva, não tumor primário.)
#   C. Extrai o transcriptoma completo DAS AMOSTRAS DE FÍGADO (streaming gzip+cut).
#   D. Converte Ensembl ID -> símbolo HUGO (org.Hs.eg.db) e agrega isoformas
#      (média por símbolo). Gera:
#        - liver_transcriptome_symbols.tsv  (transcriptoma p/ GSVA, símbolos)
#        - liver_panel_41genes.tsv          (painel ECM/MMP/TIMP, símbolos)
#
# JUSTIFICATIVA DE BAIXAR O TRANSCRIPTOMA COMPLETO:
#   (i)  Os endpoints de single-gene do hub Toil respondem HTTP 403 (S3); só os
#        arquivos físicos .gz estão acessíveis -> é preciso baixar o arquivo
#        completo e filtrar localmente.
#   (ii) GSVA exige o ranking de TODOS os genes por amostra como fundo.
#
# NOTA SOBRE OS IDs: o dataset TcgaTargetGtex_rsem_gene_tpm usa Ensembl gene IDs
#   (ex.: ENSG00000100985.16). Por isso o mapeamento para símbolo é obrigatório.
#
# Executar a partir da pasta nanotec2026/:
#   Rscript scripts/00_download_data.R
# =============================================================================

cat("============================================================\n")
cat("NANOTEC 2026 - Etapa 00: Download + extração + mapeamento\n")
cat("============================================================\n")

base     <- "."
dir_raw  <- file.path(base, "dados", "raw")
dir_proc <- file.path(base, "dados", "processed")
dir_docs <- file.path(base, "docs")
dir.create(dir_raw,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_proc, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_docs, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Pacotes (instalação automática se ausente)
# -----------------------------------------------------------------------------
pkg_cran <- c("data.table")
pkg_bioc <- c("UCSCXenaTools", "org.Hs.eg.db", "AnnotationDbi")

for (p in pkg_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (p in pkg_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages(library(UCSCXenaTools))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(data.table))
cat("[OK] Pacotes carregados (UCSCXenaTools", as.character(packageVersion("UCSCXenaTools")), ").\n\n")

# -----------------------------------------------------------------------------
# 2. Painel de genes ECM/MMP/TIMP (protocolo)
# -----------------------------------------------------------------------------
panel <- c(
  "MMP1","MMP2","MMP3","MMP7","MMP8","MMP9","MMP10","MMP11","MMP12","MMP13",
  "MMP14","MMP15","MMP16",
  "TIMP1","TIMP2","TIMP3","TIMP4",
  "COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL6A1",
  "LAMC1","LAMA1","LAMB1",
  "FN1","SPP1","THBS1",
  "ITGB1","ITGAV","ITGA5","ITGA6",
  "PLAU","PLAUR","SERPINE1","VEGFA","TGFB1","CDH1","VIM"
)
panel <- unique(panel)
cat(sprintf("[PAINEL] %d genes ECM/MMP/TIMP definidos.\n\n", length(panel)))

# -----------------------------------------------------------------------------
# 3. Download via UCSCXenaTools (skip se já existir)
# -----------------------------------------------------------------------------
datasets <- c(
  tpm   = "TcgaTargetGtex_rsem_gene_tpm",
  pheno = "TcgaTargetGTEX_phenotype.txt",
  surv  = "TCGA_survival_data"
)

download_dataset <- function(ds_name, destdir) {
  local_files <- list.files(destdir, full.names = TRUE)
  hit <- grep(gsub("[.]", "\\\\.", basename(ds_name)), local_files, value = TRUE)
  hit <- hit[file.info(hit)$size > 0]
  if (length(hit) > 0) {
    cat(sprintf("[SKIP] '%s' já presente (%s, %.1f KB).\n",
                ds_name, basename(hit[1]), file.info(hit[1])$size/1024))
    return(hit[1])
  }
  cat(sprintf("[DOWNLOAD] %s ...\n", ds_name))
  xe <- XenaGenerate(subset = XenaHostNames == "toilHub" & XenaDatasets == ds_name)
  q  <- XenaQuery(xe)
  XenaDownload(q, destdir = destdir, force = FALSE)
  Sys.sleep(1)
  local_files <- list.files(destdir, full.names = TRUE)
  hit <- grep(gsub("[.]", "\\\\.", basename(ds_name)), local_files, value = TRUE)
  hit <- hit[file.info(hit)$size > 0]
  if (length(hit) == 0) stop("Falha ao baixar: ", ds_name)
  cat(sprintf("[OK] %s baixado (%.1f MB).\n", ds_name, file.info(hit[1])$size/1e6))
  return(hit[1])
}

f_tpm   <- download_dataset(datasets["tpm"],   dir_raw)
f_pheno <- download_dataset(datasets["pheno"], dir_raw)
f_surv  <- download_dataset(datasets["surv"],  dir_raw)

# -----------------------------------------------------------------------------
# 4. Fenótipo + grupos hepáticos
# -----------------------------------------------------------------------------
cat("\n[FENOTIPO] Lendo fenótipo...\n")
read_table <- function(f) {
  if (grepl("\\.gz$", f, ignore.case = TRUE)) {
    read.delim(gzfile(f), header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               check.names = FALSE, quote = "")
  } else {
    read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               check.names = FALSE, quote = "")
  }
}
pheno <- read_table(f_pheno)
cat(sprintf("[FENOTIPO] %d amostras x %d colunas.\n", nrow(pheno), ncol(pheno)))

grp_normal   <- pheno[["_study"]] == "GTEX" & pheno[["_primary_site"]] == "Liver" &
                pheno[["_sample_type"]] == "Normal Tissue"
grp_tumor    <- pheno[["_study"]] == "TCGA" & pheno[["_primary_site"]] == "Liver" &
                pheno[["_sample_type"]] == "Primary Tumor"
grp_adjacent <- pheno[["_study"]] == "TCGA" & pheno[["_primary_site"]] == "Liver" &
                pheno[["_sample_type"]] == "Solid Tissue Normal"

pheno[["grupo"]] <- NA_character_
pheno[["grupo"]][grp_normal]   <- "Normal"
pheno[["grupo"]][grp_tumor]    <- "Tumor"
pheno[["grupo"]][grp_adjacent] <- "Adjacent"

liver <- pheno[!is.na(pheno[["grupo"]]), ]
cat(sprintf("[GRUPOS] Amostras de fígado retidas: %d\n", nrow(liver)))
print(table(liver[["grupo"]]))

annot <- data.frame(
  sample = liver[["sample"]], grupo = liver[["grupo"]],
  estudo = liver[["_study"]], tecido = liver[["_sample_type"]],
  stringsAsFactors = FALSE
)
write.table(annot, file.path(dir_proc, "samples_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Extração do transcriptoma completo das amostras de fígado (streaming)
# -----------------------------------------------------------------------------
cat("\n[EXTRACAO] Localizando amostras de fígado no arquivo TPM...\n")
con <- if (grepl("\\.gz$", f_tpm, ignore.case = TRUE)) gzfile(f_tpm, "rt") else file(f_tpm, "rt")
header_line <- readLines(con, n = 1)
close(con)
samples_tpm <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
sample_ids <- samples_tpm[-1]
cat(sprintf("[EXTRACAO] %d colunas de amostra no TPM.\n", length(sample_ids)))

sel <- which(sample_ids %in% liver[["sample"]])
col_idx <- sel + 1
cat(sprintf("[EXTRACAO] %d amostras de fígado presentes na matriz.\n", length(col_idx)))
if (length(col_idx) < 400) stop("POUCAS amostras de fígado na matriz. Abortando.")

idx_list <- paste(c(1, col_idx), collapse = ",")
out_full <- file.path(dir_proc, "liver_transcriptome_ensembl.tsv")
if (!file.exists(out_full)) {
  cmd_full <- paste0("gzip -cd \"", f_tpm, "\" | cut -f", idx_list,
                     " > \"", out_full, "\"")
  cat("[EXTRACAO] Extraindo transcriptoma (Ensembl) das amostras de fígado...\n")
  shell(cmd_full)
} else {
  cat("[SKIP] Transcriptoma Ensembl já extraído.\n")
}

# -----------------------------------------------------------------------------
# 6. Mapeamento Ensembl -> símbolo + agregação de isoformas
# -----------------------------------------------------------------------------
cat("\n[MAPA] Lendo transcriptoma (Ensembl) e convertendo para símbolo HUGO...\n")
dt <- fread(out_full, header = TRUE, sep = "\t")
setnames(dt, 1, "ensembl")
cat(sprintf("[MAPA] %d linhas (genes/versões Ensembl) x %d amostras lidas.\n",
            nrow(dt), ncol(dt) - 1))

# remove o sufixo de versão (.N) do Ensembl ID
dt[, ensembl := sub("\\..*$", "", ensembl)]

# mapeia Ensembl -> símbolo (org.Hs.eg.db, local)
map <- select(org.Hs.eg.db, keys = unique(dt[["ensembl"]]), keytype = "ENSEMBL",
              columns = "SYMBOL")
map <- map[!is.na(map[["SYMBOL"]]), c("ENSEMBL", "SYMBOL")]
cat(sprintf("[MAPA] %d Ensembl IDs mapeados para %d símbolos únicos.\n",
            length(unique(map[["ENSEMBL"]])), length(unique(map[["SYMBOL"]]))))

# junta e descarta linhas sem símbolo
dt <- merge(dt, map, by.x = "ensembl", by.y = "ENSEMBL", all.x = TRUE)
n_sem_simbolo <- sum(is.na(dt[["SYMBOL"]]))
dt <- dt[!is.na(dt[["SYMBOL"]])]
dt[, ensembl := NULL]
cat(sprintf("[MAPA] %d linhas sem símbolo descartadas; %d linhas com símbolo.\n",
            n_sem_simbolo, nrow(dt)))

# agrega isoformas/versões por símbolo (média dos TPM)
sample_cols <- setdiff(names(dt), "SYMBOL")
cat("[MAPA] Agregando linhas por símbolo (média)...\n")
dt_sym <- dt[, lapply(.SD, mean), by = SYMBOL, .SDcols = sample_cols]
cat(sprintf("[MAPA] Transcriptoma final: %d símbolos x %d amostras.\n",
            nrow(dt_sym), length(sample_cols)))

out_sym <- file.path(dir_proc, "liver_transcriptome_symbols.tsv")
fwrite(dt_sym, out_sym, sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# 7. Painel (41 genes) em símbolos
# -----------------------------------------------------------------------------
painel_dt <- dt_sym[SYMBOL %in% panel]
setcolorder(painel_dt, c("SYMBOL", sample_cols))
out_panel <- file.path(dir_proc, "liver_panel_41genes.tsv")
fwrite(painel_dt, out_panel, sep = "\t", quote = FALSE)
cat(sprintf("[PAINEL] %d genes do painel presentes após mapeamento -> %s\n",
            nrow(painel_dt), out_panel))

# -----------------------------------------------------------------------------
# 8. Manifesto + resumo
# -----------------------------------------------------------------------------
manifest <- data.frame(
  dataset = names(datasets),
  nome    = c(f_tpm, f_pheno, f_surv),
  tamanho_bytes = sapply(c(f_tpm, f_pheno, f_surv), function(f) file.info(f)$size),
  url = c(
    "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz",
    "https://toil.xenahubs.net/download/TCGA_survival_data"
  ),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(dir_docs, "manifesto_download.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

resumo <- data.frame(
  item = c("Amostras fígado (total)", "Normal (GTEx Liver)",
           "Tumor (TCGA-LIHC Primary)", "Adjacent (TCGA Solid Tissue Normal)",
           "Genes do painel (protocolo)", "Genes do painel na matriz",
           "Transcriptoma (símbolos únicos)"),
  valor = c(nrow(liver), sum(grp_normal), sum(grp_tumor), sum(grp_adjacent),
            length(panel), nrow(painel_dt), nrow(dt_sym)),
  stringsAsFactors = FALSE
)
write.table(resumo, file.path(dir_docs, "resumo_download.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n[FIM] Etapa 00 concluída.\n")
cat("  Painel        :", out_panel, "\n")
cat("  Transcriptoma :", out_sym, "\n")
cat("  Anotação      :", file.path(dir_proc, "samples_annotation.tsv"), "\n")
