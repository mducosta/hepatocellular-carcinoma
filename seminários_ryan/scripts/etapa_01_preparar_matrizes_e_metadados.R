# =============================================================================
# etapa_01_preparar_matrizes_e_metadados.R
# -----------------------------------------------------------------------------
# Seminário Ryan — Hepatocarcinoma (TCGA-LIHC) x GTEx (fígado)
#
# Objetivo:
#   1. Ler os dados de expressão (transcriptoma completo do hub Toil/Xena,
#      quando disponível; caso contrário, o painel hepático liver.tsv).
#   2. Construir a matriz de expressão em log2(TPM + 0,001), com símbolos HUGO.
#   3. Integrar fenótipos clínicos de TCGA-LIHC e GTEx (sexo, idade, estágio,
#      grau histológico, etc.).
#   4. Salvar matrizes e metadados processados para as etapas seguintes.
#
# Entrada esperada:
#   dados/raw/liver.tsv                          (painel hepático; 213 genes)
#   dados/raw/LIHC_clinicalMatrix                (clínico TCGA-LIHC — Xena)
#   dados/raw/TCGA-LIHC.clinical.tsv.gz          (clínico TCGA-LIHC — GDC)
#   dados/raw/GTEX_phenotype.gz                  (fenótipo GTEx — Xena)
#   dados/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
#   /tmp/hepato_full/TcgaTargetGtex_rsem_gene_tpm.gz  (opcional; Toil/Xena)
#
# Saídas:
#   dados/processed/metadata_amostras.tsv
#   dados/processed/metadata_analise.tsv
#   dados/processed/expressao_log2tpm_filtrada.tsv.gz
#   dados/processed/expressao_LA_pathway.tsv.gz
#   dados/processed/resumo_preparacao.txt
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# -----------------------------------------------------------------------------
# 0. Caminhos
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
REPO <- normalizePath(file.path(ROOT, ".."), mustWork = TRUE)

dir_proc <- file.path(ROOT, "dados", "processed")
dir_raw  <- file.path(ROOT, "dados", "raw")
dir.create(dir_proc, recursive = TRUE, showWarnings = FALSE)

tpm_candidates <- c(
  Sys.getenv("TCGA_GTEX_TPM", unset = ""),
  file.path(Sys.getenv("TEMP"), "hepato_full", "TcgaTargetGtex_rsem_gene_tpm.gz"),
  file.path(Sys.getenv("TMP"), "hepato_full", "TcgaTargetGtex_rsem_gene_tpm.gz")
)
tpm_candidates <- tpm_candidates[nzchar(tpm_candidates) & file.exists(tpm_candidates)]
tpm_file <- if (length(tpm_candidates) > 0) tpm_candidates[1] else ""
liver_file <- file.path(dir_raw, "liver.tsv")
if (!file.exists(liver_file)) liver_file <- file.path(REPO, "dados", "raw", "liver.tsv")
clin_file  <- file.path(dir_raw, "LIHC_clinicalMatrix")
gdc_file   <- file.path(dir_raw, "TCGA-LIHC.clinical.tsv.gz")
gtex_pheno_file <- file.path(dir_raw, "GTEX_phenotype.gz")
gtex_subj_file  <- file.path(dir_raw, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

stopifnot(file.exists(liver_file), file.exists(clin_file),
          file.exists(gtex_pheno_file), file.exists(gtex_subj_file))
if (!file.exists(gdc_file)) warning("GDC clinical não encontrado: ", gdc_file)

log_lines <- c(
  "=============================================================",
  "PREPARAÇÃO DOS DADOS — TCGA-LIHC x GTEx (fígado)",
  "=============================================================",
  paste("Data:", Sys.time()),
  paste("Raiz do seminário:", ROOT),
  paste("Repositório:", REPO)
)

# -----------------------------------------------------------------------------
# 1. Fenótipo básico a partir de liver.tsv
# -----------------------------------------------------------------------------
message("[1/5] Lendo fenótipo básico (liver.tsv) ...")
liver <- fread(liver_file, sep = "\t", data.table = FALSE, check.names = FALSE)
liver <- liver[, c("sample", "_gender", "_sample_type", "_study",
                   "TCGA_GTEX_main_category", "OS", "OS.time",
                   "PFI", "PFI.time")]
names(liver) <- c("sample", "sex_liver", "sample_type", "study",
                  "main_category", "OS", "OS_time", "PFI", "PFI_time")
liver$study <- toupper(liver$study)

# Condição biológica
liver$condition <- NA_character_
liver$condition[liver$study == "TCGA" & liver$sample_type == "Primary Tumor"] <- "Tumor"
liver$condition[liver$study == "TCGA" & liver$sample_type == "Solid Tissue Normal"] <- "Adjacent"
liver$condition[liver$study == "TCGA" & liver$sample_type == "Recurrent Tumor"] <- "Recurrent"
liver$condition[liver$study == "GTEX" & grepl("Normal", liver$sample_type, ignore.case = TRUE)] <- "Normal"

liver$group <- ifelse(liver$condition == "Tumor", "TCGA_Tumor",
               ifelse(liver$condition == "Adjacent", "TCGA_Adjacent",
               ifelse(liver$condition == "Recurrent", "TCGA_Recurrent",
               ifelse(liver$condition == "Normal", "GTEx_Normal", NA))))
liver$dataset <- ifelse(liver$study == "TCGA", "TCGA-LIHC", "GTEx-Liver")

# -----------------------------------------------------------------------------
# 2. Clínico TCGA-LIHC (Xena/TCGA)
# -----------------------------------------------------------------------------
message("[2/5] Integrando clínico TCGA-LIHC ...")
clin <- fread(clin_file, sep = "\t", data.table = FALSE, check.names = FALSE)
names(clin)[1] <- "sample"

# Harmonização de estágio AJCC (I, II, III, IV)
harmonize_stage <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("^Stage I($|A|B)", x, ignore.case = TRUE)]  <- "I"
  out[grepl("^Stage II($|A|B)", x, ignore.case = TRUE)] <- "II"
  out[grepl("^Stage III", x, ignore.case = TRUE)]       <- "III"
  out[grepl("^Stage IV", x, ignore.case = TRUE)]        <- "IV"
  out[grepl("Discrepancy", x, ignore.case = TRUE)]      <- NA_character_
  out
}

# Harmonização de grau histológico (G1–G4)
harmonize_grade <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^G1", x)] <- "G1"
  out[grepl("^G2", x)] <- "G2"
  out[grepl("^G3", x)] <- "G3"
  out[grepl("^G4", x)] <- "G4"
  out
}

# Escore de Ishak: a coluna original é texto ("0 - No Fibrosis",
# "1,2 - Portal Fibrosis", ..., "6 - Established Cirrhosis").
parse_ishak <- function(x) {
  x <- trimws(as.character(x))
  suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", x)))
}

clin2 <- data.frame(
  sample = clin$sample,
  sex_tcga = tolower(trimws(clin$gender)),
  age = suppressWarnings(as.numeric(clin$age_at_initial_pathologic_diagnosis)),
  stage_raw = clin$pathologic_stage,
  stage = harmonize_stage(clin$pathologic_stage),
  grade_raw = clin$neoplasm_histologic_grade,
  grade = harmonize_grade(clin$neoplasm_histologic_grade),
  histology = clin$histological_type,
  vital_status = clin$vital_status,
  days_to_death = suppressWarnings(as.numeric(clin$days_to_death)),
  days_to_last_followup = suppressWarnings(as.numeric(clin$days_to_last_followup)),
  fibrosis_ishak = parse_ishak(clin$fibrosis_ishak_score),
  fibrosis_raw = clin$fibrosis_ishak_score,
  afp = suppressWarnings(as.numeric(clin$fetoprotein_outcome_value)),
  viral_hepatitis = clin$viral_hepatitis_serology,
  child_pugh = clin$child_pugh_classification_grade,
  vascular_invasion = clin$vascular_tumor_cell_type,
  residual_tumor = clin$residual_tumor,
  neoplasm_status = clin$person_neoplasm_cancer_status,
  pathologic_t = clin$pathologic_T,
  pathologic_n = clin$pathologic_N,
  pathologic_m = clin$pathologic_M,
  stringsAsFactors = FALSE
)
clin2 <- clin2[!duplicated(clin2$sample), ]

# -----------------------------------------------------------------------------
# 2b. Clínico GDC (raça, etnia, malignidade prévia, T/N/M)
# -----------------------------------------------------------------------------
if (file.exists(gdc_file)) {
  gdc <- fread(gdc_file, sep = "\t", data.table = FALSE, check.names = FALSE)
  names(gdc)[1] <- "sample"
  gdc$sample <- sub("([0-9]{2})[A-Z]$", "\\1", gdc$sample)
  gdc2 <- data.frame(
    sample = gdc$sample,
    race = gdc[["race.demographic"]],
    ethnicity = gdc[["ethnicity.demographic"]],
    prior_malignancy = gdc[["prior_malignancy.diagnoses"]],
    prior_treatment = gdc[["prior_treatment.diagnoses"]],
    ajcc_t = gdc[["ajcc_pathologic_t.diagnoses"]],
    ajcc_n = gdc[["ajcc_pathologic_n.diagnoses"]],
    ajcc_m = gdc[["ajcc_pathologic_m.diagnoses"]],
    vital_status_gdc = gdc[["vital_status.demographic"]],
    age_gdc = suppressWarnings(as.numeric(gdc[["age_at_index.demographic"]])),
    stringsAsFactors = FALSE
  )
  gdc2 <- gdc2[!duplicated(gdc2$sample), ]
} else {
  gdc2 <- data.frame(sample = character(0))
}

# -----------------------------------------------------------------------------
# 3. Clínico GTEx (sexo e faixa etária)
# -----------------------------------------------------------------------------
message("[3/5] Integrando fenótipos GTEx ...")
gtex_subj <- fread(gtex_subj_file, sep = "\t", data.table = FALSE)
gtex_subj$sex_gtex <- ifelse(gtex_subj$SEX == 1, "male",
                       ifelse(gtex_subj$SEX == 2, "female", NA))
# Ponto médio da faixa etária
age_mid <- c("20-29" = 25, "30-39" = 35, "40-49" = 45,
             "50-59" = 55, "60-69" = 65, "70-79" = 75)
gtex_subj$age_gtex <- unname(age_mid[as.character(gtex_subj$AGE)])

# -----------------------------------------------------------------------------
# 4. Montagem do metadado final
# -----------------------------------------------------------------------------
message("[4/5] Montando metadados ...")
meta <- liver
meta$subject <- ifelse(meta$study == "GTEX",
                       sub("^(GTEX-[^-]+).*$", "\\1", meta$sample),
                       sub("^(TCGA-[^-]+-[^-]+).*$", "\\1", meta$sample))

meta <- merge(meta, clin2, by = "sample", all.x = TRUE, sort = FALSE)
if (nrow(gdc2) > 0) {
  meta <- merge(meta, gdc2, by = "sample", all.x = TRUE, sort = FALSE)
} else {
  meta$race <- NA_character_; meta$ethnicity <- NA_character_
  meta$prior_malignancy <- NA_character_; meta$prior_treatment <- NA_character_
  meta$ajcc_t <- NA_character_; meta$ajcc_n <- NA_character_; meta$ajcc_m <- NA_character_
  meta$vital_status_gdc <- NA_character_; meta$age_gdc <- NA_real_
}
meta <- merge(meta, gtex_subj[, c("SUBJID", "sex_gtex", "age_gtex")],
              by.x = "subject", by.y = "SUBJID", all.x = TRUE, sort = FALSE)

meta$sex <- ifelse(meta$study == "TCGA", meta$sex_tcga, meta$sex_gtex)
meta$sex[is.na(meta$sex) | meta$sex == ""] <- meta$sex_liver[is.na(meta$sex) | meta$sex == ""]
meta$sex <- tolower(trimws(meta$sex))
meta$sex <- ifelse(grepl("^m", meta$sex), "male",
             ifelse(grepl("^f", meta$sex), "female", NA))

meta$age <- ifelse(meta$study == "TCGA",
                   ifelse(is.na(meta$age), meta$age_gdc, meta$age),
                   meta$age_gtex)
meta$age_group <- ifelse(is.na(meta$age), NA,
                   ifelse(meta$age >= 65, ">=65", "<65"))
meta$age_group60 <- ifelse(is.na(meta$age), NA,
                     ifelse(meta$age >= 60, ">=60", "<60"))

# Estágio/grau agrupados (apenas TCGA tumor)
meta$stage_group <- ifelse(meta$stage %in% c("I", "II"), "I-II",
                     ifelse(meta$stage %in% c("III", "IV"), "III-IV", NA))
meta$grade_group <- ifelse(meta$grade %in% c("G1", "G2"), "G1-G2",
                     ifelse(meta$grade %in% c("G3", "G4"), "G3-G4", NA))

# -----------------------------------------------------------------------------
# 4b. Fenótipos derivados (vírus, AFP, fibrose, Child-Pugh, histologia)
# -----------------------------------------------------------------------------
meta$tss <- ifelse(meta$study == "TCGA",
                   sub("^TCGA-([^-]+)-.*$", "\\1", meta$sample), NA_character_)

# Marcadores sorológicos (a base registra a presença de marcadores, não
# necessariamente infecção ativa; interpretar como "marcador presente").
meta$hbv_marker <- grepl("Hepatitis B Surface Antigen|HBV", meta$viral_hepatitis,
                         ignore.case = TRUE)
meta$hcv_marker <- grepl("Hepatitis[[:space:]]+C|HCV", meta$viral_hepatitis,
                         ignore.case = TRUE)
meta$viral_status <- ifelse(is.na(meta$viral_hepatitis) | meta$viral_hepatitis == "",
                            "Unknown",
                     ifelse(meta$hbv_marker & meta$hcv_marker, "HBV+HCV",
                     ifelse(meta$hbv_marker, "HBV",
                     ifelse(meta$hcv_marker, "HCV", "None/Other"))))

meta$afp_log2 <- ifelse(is.na(meta$afp), NA_real_, log2(meta$afp + 1))
meta$afp_high <- ifelse(is.na(meta$afp), NA_character_,
                 ifelse(meta$afp >= 400, ">=400", "<400"))
meta$afp_group <- ifelse(is.na(meta$afp), NA_character_,
                  ifelse(meta$afp >= 400, "AFP >= 400",
                  ifelse(meta$afp >= 20, "AFP 20-399", "AFP < 20")))

meta$fibrosis_advanced <- ifelse(is.na(meta$fibrosis_ishak), NA_character_,
                          ifelse(meta$fibrosis_ishak >= 5, ">=5", "<5"))
meta$fibrosis_group <- ifelse(is.na(meta$fibrosis_ishak), NA_character_,
                       ifelse(meta$fibrosis_ishak == 0, "0",
                       ifelse(meta$fibrosis_ishak <= 4, "1-4",
                       ifelse(meta$fibrosis_ishak == 5, "5", "6"))))

meta$child_pugh_group <- ifelse(is.na(meta$child_pugh), NA_character_,
                         ifelse(meta$child_pugh == "A", "A", "B-C"))
meta$histology_group <- ifelse(is.na(meta$histology), NA_character_,
                        ifelse(grepl("Hepatocellular", meta$histology, ignore.case = TRUE),
                               "HCC", "Other"))
meta$vascular_invasion <- ifelse(is.na(meta$vascular_invasion), "None/Unknown",
                          meta$vascular_invasion)

# Raça/etnia em grupos analíticos
meta$race_group <- meta$race
meta$race_group[is.na(meta$race_group) | meta$race_group == "not reported"] <- "Unknown"
meta$ethnicity_group <- meta$ethnicity
meta$ethnicity_group[is.na(meta$ethnicity_group) | meta$ethnicity_group == "not reported"] <- "Unknown"

# Ordem das amostras e dos fatores
meta <- meta[match(liver$sample, meta$sample), ]
meta$condition <- factor(meta$condition,
                         levels = c("Tumor", "Adjacent", "Normal", "Recurrent"))
meta$group <- factor(meta$group,
                     levels = c("TCGA_Tumor", "TCGA_Adjacent", "GTEx_Normal", "TCGA_Recurrent"))
meta$study <- factor(meta$study, levels = c("TCGA", "GTEX"))
meta$sex <- factor(meta$sex, levels = c("female", "male"))
meta$stage <- factor(meta$stage, levels = c("I", "II", "III", "IV"))
meta$grade <- factor(meta$grade, levels = c("G1", "G2", "G3", "G4"))
meta$viral_status <- factor(meta$viral_status,
                            levels = c("None/Other", "HBV", "HCV", "HBV+HCV", "Unknown"))
meta$child_pugh <- factor(meta$child_pugh, levels = c("A", "B", "C"))
meta$vascular_invasion <- factor(meta$vascular_invasion,
                                 levels = c("None/Unknown", "Micro", "Macro"))
meta$afp_high <- factor(meta$afp_high, levels = c("<400", ">=400"))
meta$fibrosis_advanced <- factor(meta$fibrosis_advanced, levels = c("<5", ">=5"))
meta$histology_group <- factor(meta$histology_group, levels = c("HCC", "Other"))

# Metadado para análises principais (exclui recidiva, n = 2)
meta_analise <- meta[meta$condition != "Recurrent", ]
meta_analise$condition <- droplevels(meta_analise$condition)
meta_analise$group <- droplevels(meta_analise$group)

write.table(meta, file.path(dir_proc, "metadata_amostras.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(meta_analise, file.path(dir_proc, "metadata_analise.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# -----------------------------------------------------------------------------
# 5. Matriz de expressão
# -----------------------------------------------------------------------------
message("[5/5] Preparando matriz de expressão ...")

if (file.exists(tpm_file)) {
  message("  Usando transcriptoma completo: ", tpm_file)
  # Lê apenas a coluna de ID e as amostras hepáticas. O cabeçalho do arquivo
  # Toil tem a primeira coluna chamada "sample" (que contém IDs Ensembl).
  con_h <- gzfile(tpm_file, "rt")
  header_line <- readLines(con_h, n = 1)
  close(con_h)
  header_cols <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  idx_samples <- which(header_cols %in% meta$sample)
  tpm <- fread(tpm_file, sep = "\t", data.table = TRUE,
               select = c(1, idx_samples), showProgress = TRUE)
  id_col <- names(tpm)[1]
  setnames(tpm, id_col, "ensembl")
  tpm[, ensembl := sub("\\..*$", "", ensembl)]

  # Mapeamento Ensembl -> símbolo HUGO
  message("  Mapeando Ensembl -> símbolo HUGO ...")
  map <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = unique(tpm$ensembl),
                               keytype = "ENSEMBL",
                               columns = c("SYMBOL"))
  map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", c("ENSEMBL", "SYMBOL")]
  map <- map[!duplicated(map$ENSEMBL), ]
  names(map) <- c("ensembl", "SYMBOL")
  map <- as.data.table(map)
  tpm <- merge(tpm, map, by = "ensembl", all.x = FALSE, sort = FALSE)
  tpm[, ensembl := NULL]

  sample_cols <- setdiff(names(tpm), "SYMBOL")
  # Agrega múltiplas linhas por símbolo pela média
  message("  Agregando isoformas/duplicatas por símbolo ...")
  mat <- as.matrix(tpm[, ..sample_cols])
  rownames(mat) <- tpm$SYMBOL
  grp <- rownames(mat)
  sums <- rowsum(mat, group = grp, reorder = FALSE)
  counts <- tabulate(match(grp, unique(grp)))
  mat <- sums / counts
  rm(tpm, sums)

  # Filtro: expressão detectável (log2 TPM >= 1) em >= 20% das amostras
  keep <- rowMeans(mat >= 1, na.rm = TRUE) >= 0.20
  mat <- mat[keep, , drop = FALSE]
  mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

  # Garante que a ordem das colunas = ordem do metadado
  common <- intersect(meta_analise$sample, colnames(mat))
  mat <- mat[, common, drop = FALSE]

  fwrite(data.table(gene = rownames(mat), as.data.frame(mat)),
         file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"),
         sep = "\t", compress = "gzip")
  n_genes_full <- nrow(mat)
  n_samples_full <- ncol(mat)
} else {
  message("  AVISO: transcriptoma completo não encontrado. Usando liver.tsv.")
  n_genes_full <- NA_integer_
  n_samples_full <- NA_integer_
}

# Painel da via KEGG hsa05417 (Lipídios e Aterosclerose) — sempre disponível
la <- fread(liver_file, sep = "\t", data.table = FALSE, check.names = FALSE)
la_meta_cols <- c("sample", "samples", "detailed_category", "DFI", "DFI.time",
                  "DSS", "DSS.time", "_gender", "TCGA_GTEX_main_category",
                  "OS", "OS.time", "PFI", "PFI.time",
                  "primary disease or tissue", "_primary_site",
                  "_sample_type", "_study")
la_genes <- setdiff(names(la), la_meta_cols)
la_expr <- as.matrix(la[, la_genes])
rownames(la_expr) <- la$sample
la_expr <- t(la_expr)
la_expr <- la_expr[, intersect(meta$sample, colnames(la_expr)), drop = FALSE]
fwrite(data.table(gene = rownames(la_expr), as.data.frame(la_expr)),
       file.path(dir_proc, "expressao_LA_pathway.tsv.gz"),
       sep = "\t", compress = "gzip")

# -----------------------------------------------------------------------------
# 6. Resumo
# -----------------------------------------------------------------------------
resumo <- c(
  log_lines,
  "",
  "AMOSTRAS (liver.tsv):",
  paste(capture.output(print(table(meta$study, meta$condition, useNA = "ifany"))), collapse = "\n"),
  "",
  paste("Total de amostras:", nrow(meta)),
  paste("Amostras para análise principal (sem recidiva):", nrow(meta_analise)),
  "",
  "TRANSCRIPTOMA COMPLETO:",
  paste("Genes após filtro:", n_genes_full),
  paste("Amostras na matriz:", n_samples_full),
  "",
  "PAINEL LA (KEGG hsa05417):",
  paste("Genes:", nrow(la_expr), "| Amostras:", ncol(la_expr)),
  "",
  "CLÍNICO TCGA:",
  paste("Estágio:", paste(names(table(meta$stage)), table(meta$stage), sep = "=", collapse = "; ")),
  paste("Grau:", paste(names(table(meta$grade)), table(meta$grade), sep = "=", collapse = "; ")),
  paste("Idade (TCGA):", paste(round(summary(meta$age[meta$study == "TCGA"]), 1), collapse = " | ")),
  paste("Idade (GTEx):", paste(round(summary(meta$age[meta$study == "GTEX"]), 1), collapse = " | ")),
  paste("Vírus (TCGA):", paste(names(table(meta$viral_status)), table(meta$viral_status), sep = "=", collapse = "; ")),
  paste("AFP>=400 (TCGA):", paste(names(table(meta$afp_high)), table(meta$afp_high), sep = "=", collapse = "; ")),
  paste("Fibrose avançada (TCGA):", paste(names(table(meta$fibrosis_advanced)), table(meta$fibrosis_advanced), sep = "=", collapse = "; ")),
  "",
  "Sessão R:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(resumo, file.path(dir_proc, "resumo_preparacao.txt"))
cat(paste(resumo, collapse = "\n"), "\n")
