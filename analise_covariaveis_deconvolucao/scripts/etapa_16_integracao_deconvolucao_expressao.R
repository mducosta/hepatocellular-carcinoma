# =============================================================================
# etapa_16_integracao_deconvolucao_expressao.R
# -----------------------------------------------------------------------------
# Integra a deconvolução (xCell) à expressão diferencial.
#
# Pergunta: quanto do sinal de expressão diferencial é atribuível à composição
# celular estimada? Ajustar por composição pode:
#   (i) revelar genes cujo sinal não é explicado por mudança de proporção;
#   (ii) remover sinal biológico quando a composição é mediadora ou consequência
#        do próprio grupo (sobreajuste).
#
# Por isso, os modelos ajustados por composição são apresentados como análise
# de atenuação, não como modelo principal.
#
# Covariáveis de composição: xCell Hepatocytes, ImmuneScore e StromaScore
# (escores independentes de proporção, mas correlacionados com o grupo).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(car)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))

dir_proc <- file.path(ROOT, "dados", "processed")
dir_out  <- file.path(ROOT, "resultados", "analises_avancadas")
dir_tab  <- file.path(dir_out, "tabelas", "integracao_composicao")
dir_qc   <- file.path(dir_out, "qc")
for (d in c(dir_tab, dir_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
FDR_CUT <- 0.05
LFC_CUT <- 1

X <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))
uni <- fread(file.path(dir_out, "tabelas", "deconvolucao",
                       "deconvolucao_unificada_amostras.csv"),
             sep = ",", data.table = FALSE, na.strings = c("", "NA"))
uni <- uni[uni$sample %in% colnames(X), , drop = FALSE]
for (v in c("group", "sex", "stage_group", "grade_group", "viral_status",
            "child_pugh_group", "histology_group", "tss_collapsed",
            "afp_group2", "fibrosis_adv2")) {
  if (v %in% names(uni)) uni[[v]] <- factor(uni[[v]])
}
uni$stage_num <- as.numeric(factor(uni$stage, levels = c("I", "II", "III", "IV")))
uni$grade_num <- as.numeric(factor(uni$grade, levels = c("G1", "G2", "G3", "G4")))

# Escores de composição padronizados
comp_vars <- c("xCell_Hepatocytes", "xCell_ImmuneScore", "xCell_StromaScore")
for (v in comp_vars) {
  z <- as.numeric(scale(uni[[v]]))
  uni[[paste0("z_", v)]] <- z
}
comp_z <- paste0("z_", comp_vars)

# Diagnóstico de colinearidade
diag_vif <- list()
vars <- comp_z
if (length(vars) >= 2) {
  f <- as.formula(paste("~", paste(vars, collapse = " + ")))
  mm <- model.matrix(f, data = uni[complete.cases(uni[, vars]), ])
  vif <- tryCatch(car::vif(mm), error = function(e) NULL)
  if (!is.null(vif)) {
    diag_vif[["composicao"]] <- data.frame(
      conjunto = paste(vars, collapse = "+"),
      variavel = rownames(vif), vif = vif[, 1], stringsAsFactors = FALSE)
  }
}
if (length(diag_vif)) fwrite(do.call(rbind, diag_vif),
                             file.path(dir_tab, "diagnostico_VIF_composicao.csv"))

# Executa base vs composição
rodar_par <- function(nome, data, group_col, levels, cov_base,
                      cov_comp = comp_z, note = "") {
  message("Rodando: ", nome)
  base <- rodar_limma(X, data, group_col, levels, cov_base,
                      fdr_cut = FDR_CUT, lfc_cut = LFC_CUT)
  adj <- rodar_limma(X, data, group_col, levels, c(cov_base, cov_comp),
                     fdr_cut = FDR_CUT, lfc_cut = LFC_CUT)
  if (is.null(base) || is.null(adj)) return(NULL)

  fwrite(base, file.path(dir_tab, paste0(nome, "_base.csv")))
  fwrite(adj, file.path(dir_tab, paste0(nome, "_ajustado_composicao.csv")))

  met <- metricas_concordancia(base, adj)
  sig <- base[base$regulation != "NS", ]
  comuns <- merge(sig[, c("gene", "logFC")], adj[, c("gene", "logFC")],
                  by = "gene", suffixes = c("_base", "_adj"))
  atenuacao <- if (nrow(comuns)) {
    median(1 - abs(comuns$logFC_adj) / pmax(abs(comuns$logFC_base), 1e-9), na.rm = TRUE)
  } else NA_real_
  data.frame(
    comparacao = nome, contraste = paste(levels[2], "vs", levels[1]),
    n_base = nrow(base), n_ajustado = nrow(adj),
    degs_base = met$degs_antes, degs_ajustado = met$degs_depois,
    jaccard = met$jaccard, spearman = met$spearman_logfc,
    atenuacao_mediana_logFC = atenuacao,
    genes_base = paste(sig$gene, collapse = ";"),
    observacao = note, stringsAsFactors = FALSE
  )
}

d_tumor <- uni[uni$group == "TCGA_Tumor", , drop = FALSE]
d_adj <- uni[uni$group == "TCGA_Adjacent", , drop = FALSE]
d_gtex <- uni[uni$group == "GTEx_Normal", , drop = FALSE]
d_tcga <- uni[uni$group %in% c("TCGA_Tumor", "TCGA_Adjacent"), , drop = FALSE]
d_tg <- uni[uni$group %in% c("TCGA_Tumor", "GTEx_Normal"), , drop = FALSE]

pares <- list(
  list(nome = "tumor_vs_adjacente_composicao", data = d_tcga,
       group_col = "group", levels = c("TCGA_Adjacent", "TCGA_Tumor"),
       cov_base = c("sex", "age"),
       note = "Composição pode ser mediadora; interpretar como atenuação."),
  list(nome = "tumor_vs_gtex_composicao", data = d_tg,
       group_col = "group", levels = c("GTEx_Normal", "TCGA_Tumor"),
       cov_base = c("sex", "age"),
       note = "Grupo e composição confundidos com estudo; não causal."),
  list(nome = "estagio_composicao", data = d_tumor,
       group_col = "stage_group", levels = c("I-II", "III-IV"),
       cov_base = c("sex", "age", "grade_group"),
       note = "Estágio ajustado por grau; composição como mediadora."),
  list(nome = "estagio_clinico_composicao", data = d_tumor,
       group_col = "stage_group", levels = c("I-II", "III-IV"),
       cov_base = c("sex", "age", "grade_group", "viral_status",
                    "child_pugh_group", "histology_group"),
       note = "Estágio com ajuste clínico + composição."),
  list(nome = "grau_composicao", data = d_tumor,
       group_col = "grade_group", levels = c("G1-G2", "G3-G4"),
       cov_base = c("sex", "age", "stage_group"),
       note = "Grau ajustado por estágio; composição como mediadora."),
  list(nome = "grau_clinico_composicao", data = d_tumor,
       group_col = "grade_group", levels = c("G1-G2", "G3-G4"),
       cov_base = c("sex", "age", "stage_group", "viral_status",
                    "child_pugh_group", "histology_group"),
       note = "Grau com ajuste clínico + composição."),
  list(nome = "sexo_gtex_composicao", data = d_gtex,
       group_col = "sex", levels = c("female", "male"),
       cov_base = c("age"),
       note = "Sexo no fígado normal; testa se diferença é composicional."),
  list(nome = "sexo_tumor_composicao", data = d_tumor,
       group_col = "sex", levels = c("female", "male"),
       cov_base = c("age", "stage_group", "grade_group"),
       note = "Sexo no tumor; ajuste clínico mínimo + composição."),
  list(nome = "sexo_tumor_clinico_composicao", data = d_tumor,
       group_col = "sex", levels = c("female", "male"),
       cov_base = c("age", "stage_group", "grade_group", "viral_status",
                    "histology_group", "tss_collapsed"),
       note = "Sexo no tumor com ajuste clínico completo + composição."),
  list(nome = "sexo_adjacente_composicao", data = d_adj,
       group_col = "sex", levels = c("female", "male"),
       cov_base = c("age"),
       note = "Sexo no adjacente; n = 50, baixo poder.")
)

res <- list()
for (p in pares) {
  out <- tryCatch(do.call(rodar_par, p), error = function(e) NULL)
  if (!is.null(out)) res[[p$nome]] <- out
}
resumo <- do.call(rbind, res)
rownames(resumo) <- NULL
fwrite(resumo, file.path(dir_tab, "comparacao_base_vs_composicao.csv"))

writeLines(c(
  "INTEGRAÇÃO DECONVOLUÇÃO × EXPRESSÃO — LOG",
  "==========================================",
  paste("Data:", Sys.time()),
  paste("R:", R.version.string),
  paste("Amostras:", nrow(uni), "| Genes:", nrow(X)),
  paste("Covariáveis de composição:", paste(comp_vars, collapse = ", ")),
  "",
  paste0(resumo$comparacao, ": DEGs base=", resumo$degs_base,
         " ajustado=", resumo$degs_ajustado,
         " | Jaccard=", round(resumo$jaccard, 3),
         " | Spearman=", round(resumo$spearman, 3),
         " | atenuação mediana=", round(resumo$atenuacao_mediana_logFC, 3))
), file.path(dir_qc, "log_integracao_composicao.txt"))

cat("Concluído.\n")
print(resumo[, c("comparacao", "n_base", "n_ajustado", "degs_base",
                 "degs_ajustado", "jaccard", "spearman", "atenuacao_mediana_logFC")])
