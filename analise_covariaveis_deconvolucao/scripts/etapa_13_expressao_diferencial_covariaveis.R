# =============================================================================
# etapa_13_expressao_diferencial_covariaveis.R
# -----------------------------------------------------------------------------
# Expressão diferencial com covariáveis clínicas ampliadas.
#
# Objetivo: avaliar até onde o sinal transcriptômico permanece após ajuste por
# sexo, idade, estágio, grau, sorologia viral, Child-Pugh, histologia, AFP,
# fibrose, invasão vascular e centro TCGA (TSS).
#
# Correções importantes em relação à versão anterior:
#   - Child-Pugh ausente ("") NÃO é recodificado como B-C. A versão anterior
#     tratava 132 tumores sem Child-Pugh como B-C; aqui eles ficam NA.
#   - GTEx com dois fragmentos do mesmo doador (GTEX-WYVS) é deduplicado na
#     análise principal; a versão com 110 amostras fica como sensibilidade.
#   - "Unknown" em sorologia é uma categoria informativa, não "sem vírus".
#   - Não existem casos "None/Other" conhecidos no TCGA-LIHC; portanto, não se
#     compara "vírus positivo vs. negativo", apenas grupos sorológicos entre si
#     e "sorologia conhecida vs. desconhecida".
#
# Método: limma (lmFit + eBayes trend/robust) em log2(TPM + 0,001).
# Significância: FDR (BH) < 0,05 e |logFC| > 1.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)

dir_proc <- file.path(ROOT, "dados", "processed")
dir_out  <- file.path(ROOT, "resultados", "analises_avancadas")
dir_tab  <- file.path(dir_out, "tabelas", "expressao_diferencial_covariaveis")
dir_qc   <- file.path(dir_out, "qc")
dir_der  <- file.path(dir_out, "dados_derivados")
for (d in c(dir_tab, dir_qc, dir_der)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

FDR_CUT <- 0.05
LFC_CUT <- 1
MIN_N   <- 6

set.seed(2026)

# -----------------------------------------------------------------------------
# 1. Metadados e correções
# -----------------------------------------------------------------------------
meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
              sep = "\t", data.table = FALSE, na.strings = c("", "NA"))

# 1a. Child-Pugh: corrigir missing codificado como string vazia
meta$child_pugh <- trimws(as.character(meta$child_pugh))
meta$child_pugh_group <- ifelse(
  is.na(meta$child_pugh) | meta$child_pugh == "",
  NA_character_,
  ifelse(meta$child_pugh == "A", "A", "B-C")
)
meta$child_pugh_group <- factor(meta$child_pugh_group, levels = c("A", "B-C"))

# 1b. Fatores clínicos
meta$group <- factor(meta$group, levels = c("TCGA_Tumor", "TCGA_Adjacent", "GTEx_Normal"))
meta$condition <- factor(meta$condition, levels = c("Tumor", "Adjacent", "Normal"))
meta$sex <- factor(meta$sex, levels = c("female", "male"))
meta$sex_liver <- factor(meta$sex_liver, levels = c("Female", "Male"))
meta$stage <- factor(meta$stage, levels = c("I", "II", "III", "IV"))
meta$grade <- factor(meta$grade, levels = c("G1", "G2", "G3", "G4"))
meta$stage_group <- factor(
  ifelse(is.na(meta$stage), NA, ifelse(meta$stage %in% c("I", "II"), "I-II", "III-IV")),
  levels = c("I-II", "III-IV")
)
meta$grade_group <- factor(
  ifelse(is.na(meta$grade), NA, ifelse(meta$grade %in% c("G1", "G2"), "G1-G2", "G3-G4")),
  levels = c("G1-G2", "G3-G4")
)
meta$stage_num <- as.numeric(meta$stage)
meta$grade_num <- as.numeric(meta$grade)

# 1c. Sorologia: manter "Unknown" como categoria explícita
meta$viral_status <- trimws(as.character(meta$viral_status))
meta$viral_status <- ifelse(meta$viral_status == "", "Unknown", meta$viral_status)
meta$viral_status <- factor(meta$viral_status,
                            levels = c("HBV", "HCV", "HBV+HCV", "Unknown"))
meta$viral_known <- factor(
  ifelse(meta$viral_status == "Unknown", NA, as.character(meta$viral_status)),
  levels = c("HBV", "HCV", "HBV+HCV")
)
meta$viral_any <- factor(
  ifelse(meta$viral_status == "Unknown", NA,
         ifelse(meta$viral_status == "HBV+HCV", "HBV+HCV", "HBV ou HCV")),
  levels = c("HBV ou HCV", "HBV+HCV")
)

# 1d. Outras variáveis
meta$histology_group <- factor(
  ifelse(is.na(meta$histology) | meta$histology == "", NA,
         ifelse(grepl("Hepatocellular Carcinoma", meta$histology, fixed = TRUE), "HCC", "Other")),
  levels = c("HCC", "Other")
)
meta$afp_group2 <- factor(
  ifelse(is.na(meta$afp), NA,
         ifelse(meta$afp >= 400, "AFP >= 400", "AFP < 400")),
  levels = c("AFP < 400", "AFP >= 400")
)
meta$fibrosis_adv2 <- factor(
  ifelse(is.na(meta$fibrosis_ishak), NA,
         ifelse(meta$fibrosis_ishak >= 5, "Ishak >= 5", "Ishak < 5")),
  levels = c("Ishak < 5", "Ishak >= 5")
)
meta$vascular_invasion2 <- factor(
  ifelse(is.na(meta$vascular_invasion) | meta$vascular_invasion == "", NA,
         as.character(meta$vascular_invasion)),
  levels = c("Micro", "Macro")
)
meta$age <- as.numeric(meta$age)

# 1e. TSS: agrupar centros pequenos para evitar soluções singulares
tss_counts <- table(meta$tss[meta$group != "GTEx_Normal"])
tss_keep <- names(tss_counts)[tss_counts >= 10]
meta$tss_collapsed <- factor(
  ifelse(meta$group == "GTEx_Normal", "GTEx",
         ifelse(meta$tss %in% tss_keep, as.character(meta$tss), "Outros")),
  levels = c(sort(tss_keep), "Outros", "GTEx")
)

# 1f. Deduplicação GTEx: um doador com dois fragmentos de fígado
meta$duplicata_gtex <- duplicated(meta$subject) & meta$group == "GTEx_Normal"
meta_prim <- meta[!meta$duplicata_gtex, , drop = FALSE]

# Fatores com níveis não usados após subset
clean_factors <- function(d) {
  for (v in names(d)) if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  d
}
meta_prim <- clean_factors(meta_prim)

d_tumor  <- meta_prim[meta_prim$group == "TCGA_Tumor", , drop = FALSE]
d_adj    <- meta_prim[meta_prim$group == "TCGA_Adjacent", , drop = FALSE]
d_gtex   <- meta_prim[meta_prim$group == "GTEx_Normal", , drop = FALSE]
d_tcga   <- meta_prim[meta_prim$group %in% c("TCGA_Tumor", "TCGA_Adjacent"), , drop = FALSE]
d_tcga_tg <- meta_prim[meta_prim$group %in% c("TCGA_Tumor", "GTEx_Normal"), , drop = FALSE]
# Subconjuntos com casos completos para comparar restrição amostral vs. ajuste
d_stage_clin <- d_tumor[complete.cases(d_tumor[, c("stage_group", "sex", "age",
  "grade_group", "viral_status", "child_pugh_group", "histology_group")]), , drop = FALSE]
d_grade_clin <- d_tumor[complete.cases(d_tumor[, c("grade_group", "sex", "age",
  "stage_group", "viral_status", "child_pugh_group", "histology_group")]), , drop = FALSE]
for (dd in list(d_tumor, d_adj, d_gtex, d_tcga, d_tcga_tg, d_stage_clin, d_grade_clin)) dd <- clean_factors(dd)

# -----------------------------------------------------------------------------
# 2. Matriz de expressão
# -----------------------------------------------------------------------------
expr_file <- file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz")
expr_dt <- fread(expr_file, sep = "\t", data.table = FALSE)
rownames(expr_dt) <- expr_dt$gene
expr <- as.matrix(expr_dt[, setdiff(names(expr_dt), "gene"), drop = FALSE])
storage.mode(expr) <- "numeric"
rm(expr_dt)
stopifnot(all(meta_prim$sample %in% colnames(expr)))

# -----------------------------------------------------------------------------
# 3. Função de modelo
# -----------------------------------------------------------------------------
empty_summary <- function(name, group_col, levels, covariates, n1, n2, status, note) {
  data.frame(
    modelo = name,
    variavel = group_col,
    contraste = paste(levels[2], "vs", levels[1]),
    n_referencia = n1,
    n_comparacao = n2,
    n_total = n1 + n2,
    n_covariaveis = length(covariates),
    df_residual = NA_integer_,
    posto_completo = NA,
    genes_testados = NA_integer_,
    up = NA_integer_,
    down = NA_integer_,
    total_degs = NA_integer_,
    menor_fdr = NA_real_,
    status = status,
    observacao = note,
    stringsAsFactors = FALSE
  )
}

run_model <- function(name, group_col, levels, covariates = character(),
                      data = meta_prim, block_col = NULL, min_n = MIN_N,
                      note = "") {
  is_numeric_group <- is.numeric(data[[group_col]])
  if (is_numeric_group) {
    d <- data[!is.na(data[[group_col]]), , drop = FALSE]
  } else {
    d <- data[data[[group_col]] %in% levels, , drop = FALSE]
    d[[group_col]] <- factor(as.character(d[[group_col]]), levels = levels)
  }

  # Remove covariáveis constantes no subset
  if (length(covariates) > 0) {
    keep_cov <- vapply(covariates, function(v) {
      x <- d[[v]]
      length(unique(x[!is.na(x)])) > 1
    }, logical(1))
    covariates <- covariates[keep_cov]
  }

  keep <- complete.cases(d[, c(group_col, covariates, block_col), drop = FALSE])
  d <- d[keep, , drop = FALSE]
  d <- clean_factors(d)

  if (is_numeric_group) {
    n1 <- sum(d[[group_col]] == min(d[[group_col]], na.rm = TRUE))
    n2 <- sum(d[[group_col]] == max(d[[group_col]], na.rm = TRUE))
  } else {
    n1 <- sum(d[[group_col]] == levels[1])
    n2 <- sum(d[[group_col]] == levels[2])
  }
  if (nrow(d) < min_n || n1 < 2 || n2 < 2) {
    return(list(summary = empty_summary(name, group_col, levels, covariates,
                                        n1, n2, "nao_executado",
                                        "n insuficiente ou grupo com <2 amostras"),
                table = NULL))
  }

  X <- expr[, d$sample, drop = FALSE]

  if (!is.null(block_col)) {
    # Análise pareada por efeitos fixos de paciente. É equivalente a um teste
    # pareado e evita a lentidão de duplicateCorrelation em 12 mil genes.
    d[[block_col]] <- factor(d[[block_col]])
    fml <- as.formula(paste("~ 0 +", group_col, "+", block_col))
    design <- model.matrix(fml, data = d)
    colnames(design) <- sub(group_col, "group", colnames(design), fixed = TRUE)
    if (!all(paste0("group", levels) %in% colnames(design))) {
      return(list(summary = empty_summary(name, group_col, levels, covariates,
                                          n1, n2, "falhou",
                                          "colunas de grupo ausentes no design"),
                  table = NULL))
    }
    rank_ok <- qr(design)$rank == ncol(design)
    if (!rank_ok) {
      return(list(summary = empty_summary(name, group_col, levels, covariates,
                                          n1, n2, "nao_executado",
                                          "design pareado com posto incompleto"),
                  table = NULL))
    }
    fit <- lmFit(X, design)
    cvec <- setNames(rep(0, ncol(design)), colnames(design))
    cvec[paste0("group", levels[1])] <- -1
    cvec[paste0("group", levels[2])] <-  1
    fit2 <- contrasts.fit(fit, cvec)
    df_res <- fit$df.residual[1]
    corr_block <- NA_real_
  } else {
    fml <- as.formula(paste("~ 0 +", group_col,
                            if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""))
    design <- model.matrix(fml, data = d)
    colnames(design) <- sub(group_col, "group", colnames(design), fixed = TRUE)
    rank_ok <- qr(design)$rank == ncol(design)
    if (!rank_ok) {
      return(list(summary = empty_summary(name, group_col, levels, covariates,
                                          n1, n2, "nao_executado",
                                          "design com posto incompleto (colinearidade)"),
                  table = NULL))
    }
    fit <- lmFit(X, design)
    cvec <- setNames(rep(0, ncol(design)), colnames(design))
    if (is_numeric_group) {
      coef_num <- c(paste0("group", group_col), "group")
      coef_num <- coef_num[coef_num %in% colnames(design)][1]
      if (is.na(coef_num)) {
        return(list(summary = empty_summary(name, group_col, levels, covariates,
                                            n1, n2, "falhou",
                                            "coeficiente numérico ausente no design"),
                    table = NULL))
      }
      cvec[coef_num] <- 1
    } else {
      if (!all(paste0("group", levels) %in% colnames(design))) {
        return(list(summary = empty_summary(name, group_col, levels, covariates,
                                            n1, n2, "falhou",
                                            "colunas de grupo ausentes no design"),
                    table = NULL))
      }
      cvec[paste0("group", levels[1])] <- -1
      cvec[paste0("group", levels[2])] <-  1
    }
    fit2 <- contrasts.fit(fit, cvec)
    df_res <- fit$df.residual[1]
    corr_block <- NA_real_
  }

  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt$modelo <- name
  tt$contraste <- paste(levels[2], "vs", levels[1])
  tt$regulation <- ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC >  LFC_CUT, "Up",
                   ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC < -LFC_CUT, "Down", "NS"))
  tt <- tt[, c("gene", "modelo", "contraste", "logFC", "AveExpr", "t",
               "P.Value", "adj.P.Val", "B", "regulation")]
  tt <- tt[order(tt$P.Value), ]

  sm <- data.frame(
    modelo = name,
    variavel = group_col,
    contraste = paste(levels[2], "vs", levels[1]),
    n_referencia = n1,
    n_comparacao = n2,
    n_total = nrow(d),
    n_covariaveis = length(covariates),
    df_residual = df_res,
    posto_completo = rank_ok,
    genes_testados = nrow(tt),
    up = sum(tt$regulation == "Up"),
    down = sum(tt$regulation == "Down"),
    total_degs = sum(tt$regulation != "NS"),
    menor_fdr = min(tt$adj.P.Val),
    status = "ok",
    observacao = paste0(note,
                        ifelse(is.na(corr_block), "",
                               sprintf(" | correlação intraclasse = %.3f", corr_block))),
    stringsAsFactors = FALSE
  )
  list(summary = sm, table = tt)
}

# -----------------------------------------------------------------------------
# 4. Lista de modelos
# -----------------------------------------------------------------------------
models <- list(
  # --- Tumor primário vs. tecido adjacente (TCGA) ---
  list(name = "tumor_vs_adjacente_base", group_col = "group",
       levels = c("TCGA_Adjacent", "TCGA_Tumor"), data = d_tcga,
       covariates = c("sex", "age"),
       note = "Modelo basal não pareado."),
  list(name = "tumor_vs_adjacente_estagio_grau", group_col = "group",
       levels = c("TCGA_Adjacent", "TCGA_Tumor"), data = d_tcga,
       covariates = c("sex", "age", "stage_group", "grade_group"),
       note = "Ajuste por estágio e grau; casos completos."),
  list(name = "tumor_vs_adjacente_clinico_completo", group_col = "group",
       levels = c("TCGA_Adjacent", "TCGA_Tumor"), data = d_tcga,
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "viral_status", "child_pugh_group", "histology_group"),
       note = "Ajuste clínico completo; perda de adjacentes por missingness."),
  list(name = "tumor_vs_adjacente_clinico_tss", group_col = "group",
       levels = c("TCGA_Adjacent", "TCGA_Tumor"), data = d_tcga,
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "viral_status", "child_pugh_group", "histology_group",
                      "tss_collapsed"),
       note = "Ajuste clínico + centro TCGA; verificar confundimento TSS-grupo."),
  list(name = "tumor_vs_adjacente_pareado", group_col = "group",
       levels = c("TCGA_Adjacent", "TCGA_Tumor"), data = d_tcga,
       covariates = character(), block_col = "subject",
       note = "Análise pareada por paciente (efeitos fixos de paciente)."),

  # --- Tumor primário vs. GTEx ---
  list(name = "tumor_vs_gtex_base", group_col = "group",
       levels = c("GTEx_Normal", "TCGA_Tumor"), data = d_tcga_tg,
       covariates = c("sex", "age"),
       note = "Comparação estruturalmente confundida com estudo/lote."),
  list(name = "tumor_vs_gtex_sem_ajuste", group_col = "group",
       levels = c("GTEx_Normal", "TCGA_Tumor"), data = d_tcga_tg,
       covariates = character(),
       note = "Sem ajuste; usado como referência de sensibilidade."),

  # --- Estágio (dentro do tumor) ---
  list(name = "estagio_base", data = d_tumor, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age"),
       note = "Estágio III-IV vs. I-II; ajuste mínimo."),
  list(name = "estagio_clinico", data = d_tumor, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age", "grade_group", "viral_status",
                      "child_pugh_group", "histology_group"),
       note = "Estágio ajustado por grau, sorologia, Child-Pugh e histologia."),
  list(name = "estagio_clinico_tss", data = d_tumor, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age", "grade_group", "viral_status",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "Estágio + centro TCGA; TSS pode absorver parte do efeito de estágio."),
  list(name = "estagio_tendencia", data = d_tumor, group_col = "stage_num",
       levels = c(1, 4),
       covariates = c("sex", "age", "grade_num", "viral_status",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "Tendência ordinal (estágio 1-4 como variável contínua)."),

  # --- Grau (dentro do tumor) ---
  list(name = "grau_base", data = d_tumor, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age"),
       note = "Grau G3-G4 vs. G1-G2; ajuste mínimo."),
  list(name = "grau_clinico", data = d_tumor, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age", "stage_group", "viral_status",
                      "child_pugh_group", "histology_group"),
       note = "Grau ajustado por estágio, sorologia, Child-Pugh e histologia."),
  list(name = "grau_clinico_tss", data = d_tumor, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age", "stage_group", "viral_status",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "Grau + centro TCGA."),
  list(name = "grau_tendencia", data = d_tumor, group_col = "grade_num",
       levels = c(1, 4),
       covariates = c("sex", "age", "stage_num", "viral_status",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "Tendência ordinal (grau 1-4 como variável contínua)."),

  # --- Sorologia (somente entre casos conhecidos) ---
  list(name = "sorologia_hbv_vs_hcv", data = d_tumor, group_col = "viral_known",
       levels = c("HCV", "HBV"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "HBV vs. HCV; não há grupo negativo conhecido."),
  list(name = "sorologia_hbvhcv_vs_hbv", data = d_tumor, group_col = "viral_known",
       levels = c("HBV", "HBV+HCV"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "HBV+HCV vs. HBV; sorologia, não infecção ativa."),
  list(name = "sorologia_hbvhcv_vs_hcv", data = d_tumor, group_col = "viral_known",
       levels = c("HCV", "HBV+HCV"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "HBV+HCV vs. HCV; n pequeno no grupo HCV."),
  list(name = "sorologia_conhecida_vs_unknown", data = d_tumor, group_col = "viral_status",
       levels = c("Unknown", "HBV"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "child_pugh_group", "histology_group", "tss_collapsed"),
       note = "HBV conhecido vs. sorologia desconhecida; missingness informativa."),

  # --- Child-Pugh corrigido ---
  list(name = "child_pugh_a_vs_bc", data = d_tumor, group_col = "child_pugh_group",
       levels = c("A", "B-C"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "viral_status", "histology_group", "tss_collapsed"),
       note = "Child-Pugh B-C (n pequeno) vs. A; missing não é B-C."),

  # --- Histologia ---
  list(name = "histologia_hcc_vs_outros", data = d_tumor, group_col = "histology_group",
       levels = c("HCC", "Other"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "viral_status", "child_pugh_group", "tss_collapsed"),
       note = "HCC clássico vs. misto/fibrolamelar; n de 'Other' = 10."),

  # --- Sexo ---
  list(name = "sexo_gtex", data = d_gtex, group_col = "sex",
       levels = c("female", "male"),
       covariates = c("age"),
       note = "GTEx fígado normal; controle positivo XIST/genes do Y."),
  list(name = "sexo_tumor", data = d_tumor, group_col = "sex",
       levels = c("female", "male"),
       covariates = c("age", "stage_group", "grade_group", "viral_status",
                      "histology_group", "tss_collapsed"),
       note = "TCGA tumor; ajuste clínico e centro."),
  list(name = "sexo_adjacente", data = d_adj, group_col = "sex",
       levels = c("female", "male"),
       covariates = c("age"),
       note = "TCGA adjacente; n = 50, baixo poder."),

  # --- Outras covariáveis clínicas ---
  list(name = "afp_400", data = d_tumor, group_col = "afp_group2",
       levels = c("AFP < 400", "AFP >= 400"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "viral_status", "tss_collapsed"),
       note = "AFP >= 400 vs. < 400; missingness 25%."),
  list(name = "fibrose_avancada", data = d_tumor, group_col = "fibrosis_adv2",
       levels = c("Ishak < 5", "Ishak >= 5"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "viral_status", "tss_collapsed"),
       note = "Ishak >= 5 vs. < 5; missingness 43%."),
  list(name = "invasao_vascular_macro_vs_micro", data = d_tumor, group_col = "vascular_invasion2",
       levels = c("Micro", "Macro"),
       covariates = c("sex", "age", "stage_group", "grade_group",
                      "viral_status", "tss_collapsed"),
       note = "Macro vs. micro; missingness 71%; não há 'sem invasão' confiável."),
  # --- Sensibilidade: restrição amostral vs. ajuste ---
  list(name = "estagio_base_mesmos_casos", data = d_stage_clin, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age"),
       note = "Mesmos casos completos do modelo clínico, mas sem ajuste clínico."),
  list(name = "estagio_sem_grau", data = d_stage_clin, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age", "viral_status", "child_pugh_group", "histology_group"),
       note = "Estágio sem ajuste por grau; avalia absorção por grau."),
  list(name = "estagio_so_grau", data = d_stage_clin, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age", "grade_group"),
       note = "Estágio ajustado apenas por grau."),
  list(name = "estagio_so_tss", data = d_stage_clin, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age", "tss_collapsed"),
       note = "Estágio ajustado apenas por centro TCGA."),
  list(name = "grau_base_mesmos_casos", data = d_grade_clin, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age"),
       note = "Mesmos casos completos do modelo clínico de grau, sem ajuste."),
  list(name = "grau_sem_estagio", data = d_grade_clin, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age", "viral_status", "child_pugh_group", "histology_group"),
       note = "Grau sem ajuste por estágio."),
  list(name = "grau_so_estagio", data = d_grade_clin, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age", "stage_group"),
       note = "Grau ajustado apenas por estágio."),
  list(name = "child_pugh_base", data = d_tumor, group_col = "child_pugh_group",
       levels = c("A", "B-C"),
       covariates = c("sex", "age"),
       note = "Child-Pugh B-C vs. A sem ajuste clínico; n de B-C = 22."),
  list(name = "histologia_base", data = d_tumor, group_col = "histology_group",
       levels = c("HCC", "Other"),
       covariates = c("sex", "age"),
       note = "Histologia sem ajuste; n de Other = 10, descritivo."),
  list(name = "histologia_estagio_grau", data = d_tumor, group_col = "histology_group",
       levels = c("HCC", "Other"),
       covariates = c("sex", "age", "stage_group", "grade_group"),
       note = "Histologia ajustada por estágio e grau; n de Other pode cair."),
  list(name = "sorologia_hbv_vs_hcv_base", data = d_tumor, group_col = "viral_known",
       levels = c("HCV", "HBV"),
       covariates = c("sex", "age"),
       note = "HBV vs. HCV sem ajuste clínico; n pequeno."),
  list(name = "sorologia_hbvhcv_vs_hbv_base", data = d_tumor, group_col = "viral_known",
       levels = c("HBV", "HBV+HCV"),
       covariates = c("sex", "age"),
       note = "HBV+HCV vs. HBV sem ajuste clínico."),
  list(name = "afp_base", data = d_tumor, group_col = "afp_group2",
       levels = c("AFP < 400", "AFP >= 400"),
       covariates = c("sex", "age"),
       note = "AFP sem ajuste clínico."),
  list(name = "fibrose_base", data = d_tumor, group_col = "fibrosis_adv2",
       levels = c("Ishak < 5", "Ishak >= 5"),
       covariates = c("sex", "age"),
       note = "Fibrose sem ajuste clínico."),

  list(name = "estagio_so_tss_full", data = d_tumor, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age", "tss_collapsed"),
       note = "Estágio na coorte completa, ajustado apenas por centro; teste de confundimento TSS."),
  list(name = "estagio_sem_cc_ed", data = d_tumor[!d_tumor$tss %in% c("CC", "ED"), ],
       group_col = "stage_group", levels = c("I-II", "III-IV"),
       covariates = c("sex", "age"),
       note = "Estágio excluindo centros CC e ED, onde Child-Pugh está ausente."),
  list(name = "estagio_sem_estagio_iv", data = d_tumor[!d_tumor$stage %in% "IV", ],
       group_col = "stage_group", levels = c("I-II", "III-IV"),
       covariates = c("sex", "age"),
       note = "Estágio III vs. I-II, excluindo estágio IV (n = 5)."),
  list(name = "grau_so_tss_full", data = d_tumor, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age", "tss_collapsed"),
       note = "Grau na coorte completa, ajustado apenas por centro."),
  list(name = "grau_sem_DD", data = d_tumor[!d_tumor$tss %in% "DD", ],
       group_col = "grade_group", levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age"),
       note = "Grau excluindo o centro DD (39% dos tumores)."),
  list(name = "grau_base_sem_estagio", data = d_tumor, group_col = "grade_group",
       levels = c("G1-G2", "G3-G4"),
       covariates = c("sex", "age", "stage_group"),
       note = "Grau na coorte completa, ajustado por estágio."),

  list(name = "estagio_x_grau_interacao", data = d_tumor, group_col = "stage_group",
       levels = c("I-II", "III-IV"),
       covariates = c("sex", "age", "grade_group", "viral_status",
                      "child_pugh_group", "histology_group",
                      "stage_group:grade_group"),
       note = "Interação estágio x grau; exploratória, poder limitado.")
)

# -----------------------------------------------------------------------------
# 5. Executar modelos
# -----------------------------------------------------------------------------
res <- list()
for (m in models) {
  message("Rodando: ", m$name)
  out <- tryCatch(
    run_model(name = m$name, group_col = m$group_col, levels = m$levels,
              covariates = m$covariates, block_col = m$block_col,
              data = m$data, note = m$note),
    error = function(e) list(summary = empty_summary(
      m$name, m$group_col, m$levels, m$covariates, NA, NA, "erro", conditionMessage(e)),
      table = NULL)
  )
  res[[m$name]] <- out
  if (!is.null(out$table)) {
    fwrite(out$table, file.path(dir_tab, paste0("modelo_", m$name, "_completo.csv")))
    sig <- out$table[out$table$regulation != "NS", ]
    fwrite(sig, file.path(dir_tab, paste0("modelo_", m$name, "_significativos.csv")))
    top <- rbind(
      head(out$table[out$table$regulation == "Up", ][order(-out$table$logFC[out$table$regulation == "Up"]), ], 20),
      head(out$table[out$table$regulation == "Down", ][order(out$table$logFC[out$table$regulation == "Down"]), ], 20)
    )
    fwrite(top, file.path(dir_tab, paste0("modelo_", m$name, "_top20.csv")))
  }
}

resumo <- do.call(rbind, lapply(res, `[[`, "summary"))
rownames(resumo) <- NULL
fwrite(resumo, file.path(dir_tab, "resumo_modelos_expressao_diferencial_covariaveis.csv"))

nao_exec <- resumo[resumo$status != "ok", ]
fwrite(nao_exec, file.path(dir_tab, "modelos_nao_executados.csv"))

# -----------------------------------------------------------------------------
# 6. Comparação entre modelos (Jaccard e Spearman)
# -----------------------------------------------------------------------------
comparar <- function(a, b) {
  ta <- res[[a]]$table; tb <- res[[b]]$table
  if (is.null(ta) || is.null(tb)) return(NULL)
  genes_comuns <- intersect(ta$gene, tb$gene)
  sa <- ta$gene[ta$regulation != "NS"]
  sb <- tb$gene[tb$regulation != "NS"]
  jac <- if (length(union(sa, sb)) == 0) NA_real_ else length(intersect(sa, sb)) / length(union(sa, sb))
  rho <- suppressWarnings(cor(ta$logFC[match(genes_comuns, ta$gene)],
                              tb$logFC[match(genes_comuns, tb$gene)],
                              method = "spearman", use = "complete.obs"))
  data.frame(modelo_a = a, modelo_b = b,
             degs_a = length(sa), degs_b = length(sb),
             degs_comuns = length(intersect(sa, sb)),
             jaccard = jac, spearman_logfc = rho,
             n_genes_comuns = length(genes_comuns),
             stringsAsFactors = FALSE)
}

pares <- list(
  c("tumor_vs_adjacente_base", "tumor_vs_adjacente_estagio_grau"),
  c("tumor_vs_adjacente_base", "tumor_vs_adjacente_clinico_completo"),
  c("tumor_vs_adjacente_base", "tumor_vs_adjacente_clinico_tss"),
  c("tumor_vs_adjacente_base", "tumor_vs_adjacente_pareado"),
  c("tumor_vs_adjacente_pareado", "tumor_vs_adjacente_clinico_completo"),
  c("tumor_vs_gtex_sem_ajuste", "tumor_vs_gtex_base"),
  c("estagio_base", "estagio_clinico"),
  c("estagio_clinico", "estagio_clinico_tss"),
  c("estagio_base", "estagio_clinico_tss"),
  c("grau_base", "grau_clinico"),
  c("grau_clinico", "grau_clinico_tss"),
  c("grau_base", "grau_clinico_tss"),
  c("estagio_base_mesmos_casos", "estagio_clinico"),
  c("estagio_base_mesmos_casos", "estagio_sem_grau"),
  c("estagio_base_mesmos_casos", "estagio_so_grau"),
  c("estagio_base_mesmos_casos", "estagio_so_tss"),
  c("grau_base_mesmos_casos", "grau_clinico"),
  c("grau_base_mesmos_casos", "grau_sem_estagio"),
  c("grau_base_mesmos_casos", "grau_so_estagio"),
  c("child_pugh_base", "child_pugh_a_vs_bc"),
  c("histologia_base", "histologia_estagio_grau"),
  c("histologia_estagio_grau", "histologia_hcc_vs_outros"),
  c("sorologia_hbv_vs_hcv_base", "sorologia_hbv_vs_hcv"),
  c("sorologia_hbvhcv_vs_hbv_base", "sorologia_hbvhcv_vs_hbv"),
  c("afp_base", "afp_400"),
  c("fibrose_base", "fibrose_avancada"),
  c("estagio_base", "estagio_so_tss_full"),
  c("estagio_base", "estagio_sem_cc_ed"),
  c("estagio_base", "estagio_sem_estagio_iv"),
  c("grau_base", "grau_so_tss_full"),
  c("grau_base", "grau_sem_DD"),
  c("grau_base", "grau_base_sem_estagio"),
  c("sexo_gtex", "sexo_tumor"),
  c("sexo_gtex", "sexo_adjacente")
)
comp <- do.call(rbind, lapply(pares, function(p) comparar(p[1], p[2])))
if (!is.null(comp)) fwrite(comp, file.path(dir_tab, "comparacao_modelos_jaccard_spearman.csv"))

# -----------------------------------------------------------------------------
# 7. Metadados usados e trilha de auditoria
# -----------------------------------------------------------------------------
fwrite(meta_prim, file.path(dir_der, "metadados_analise_primaria_corrigidos.tsv"), sep = "\t")
fwrite(meta[meta$duplicata_gtex, ], file.path(dir_der, "metadados_gtex_duplicata_excluida.tsv"), sep = "\t")

log_lines <- c(
  "EXPRESSÃO DIFERENCIAL COM COVARIÁVEIS — LOG",
  "===========================================",
  paste("Data:", Sys.time()),
  paste("R:", R.version.string),
  paste("Amostras na análise primária:", nrow(meta_prim)),
  paste("GTEx duplicado excluído:", sum(meta$duplicata_gtex)),
  paste("Genes:", nrow(expr)),
  paste("Critério: FDR <", FDR_CUT, "e |logFC| >", LFC_CUT),
  "",
  "Modelos executados:",
  paste0("  ", resumo$modelo, ": ", resumo$status, " | n = ", resumo$n_total,
         " | DEGs = ", resumo$total_degs, " | ", resumo$observacao)
)
writeLines(log_lines, file.path(dir_qc, "log_expressao_diferencial_covariaveis.txt"))

cat("Concluído. Modelos ok:", sum(resumo$status == "ok"), "de", nrow(resumo), "\n")
print(resumo[, c("modelo", "n_referencia", "n_comparacao", "n_total",
                 "up", "down", "total_degs", "status")])
