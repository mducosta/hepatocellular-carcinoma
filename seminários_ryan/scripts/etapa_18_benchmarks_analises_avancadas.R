# =============================================================================
# etapa_18_benchmarks_analises_avancadas.R
# -----------------------------------------------------------------------------
# Benchmark computacional e consolidação de métricas das análises avançadas.
# Mede tempo de leitura, limma, ComBat, EPIC, quanTIseq e xCell; registra uso de
# memória; consolida DEGs, Jaccard, Spearman e atenuação.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sva)
  library(EPIC)
  library(quantiseqr)
  library(xCell)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))

dir_proc <- file.path(ROOT, "dados", "processed")
dir_out  <- file.path(ROOT, "resultados", "analises_avancadas")
dir_aud  <- file.path(ROOT, "auditorias")
dir.create(dir_aud, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
bench <- list()
add_bench <- function(etapa, tempo, detalhe = "") {
  bench[[length(bench) + 1]] <<- data.frame(
    etapa = etapa, tempo_segundos = as.numeric(tempo[["elapsed"]]),
    detalhe = detalhe, stringsAsFactors = FALSE)
}

gc(reset = TRUE)
t0 <- Sys.time()

t <- system.time({
  prep <- preparar_metadados(file.path(dir_proc, "metadata_analise.tsv"), dedup_gtex = TRUE)
  X <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))
})
add_bench("Leitura de metadados + matriz", t,
          paste(nrow(X), "genes ×", ncol(X), "amostras"))
meta <- prep$meta
d_tcga <- prep$tcga

# -----------------------------------------------------------------------------
# limma
# -----------------------------------------------------------------------------
t <- system.time({
  r1 <- rodar_limma(X, d_tcga, "group", c("TCGA_Adjacent", "TCGA_Tumor"),
                    c("sex", "age"))
})
add_bench("limma: tumor vs adjacente basal", t,
          paste(nrow(r1), "genes;", sum(r1$regulation != "NS"), "DEGs"))

t <- system.time({
  r2 <- rodar_limma(X, prep$tumor, "stage_group", c("I-II", "III-IV"),
                    c("sex", "age", "grade_group", "viral_status",
                      "child_pugh_group", "histology_group"))
})
add_bench("limma: estágio clínico", t,
          paste(nrow(r2), "genes;", sum(r2$regulation != "NS"), "DEGs"))

# -----------------------------------------------------------------------------
# ComBat
# -----------------------------------------------------------------------------
d_cb <- prep$tcga_gtex[complete.cases(prep$tcga_gtex[, c("study", "sex", "age")]), ]
X_cb <- X[, d_cb$sample]
t <- system.time({
  invisible(sva::ComBat(dat = X_cb, batch = factor(d_cb$study),
                        mod = NULL, par.prior = TRUE))
})
add_bench("ComBat: TCGA vs GTEx sem covariável", t,
          paste(ncol(X_cb), "amostras;", nrow(X_cb), "genes"))

d_tss <- prep$tumor[complete.cases(prep$tumor[, c("tss_collapsed", "grade_group",
                                                   "stage_group", "sex", "age")]), ]
X_tss <- X[, d_tss$sample]
mod_tss <- model.matrix(~ grade_group + stage_group + sex + age, data = d_tss)
t <- system.time({
  invisible(sva::ComBat(dat = X_tss, batch = droplevels(factor(d_tss$tss_collapsed)),
                        mod = mod_tss, par.prior = TRUE))
})
add_bench("ComBat: TSS no tumor (grau + estágio)", t,
          paste(ncol(X_tss), "amostras;", nlevels(droplevels(factor(d_tss$tss_collapsed))), "lotes"))

# -----------------------------------------------------------------------------
# Deconvolução
# -----------------------------------------------------------------------------
X_meta <- X[, meta$sample, drop = FALSE]
tpm <- log2tpm_para_tpm(X_meta)
t <- system.time({
  invisible(EPIC::EPIC(bulk = tpm, mRNA_cell = NULL, scaleExprs = TRUE,
                       withOtherCells = TRUE))
})
add_bench("EPIC", t, paste(ncol(tpm), "amostras"))

t <- system.time({
  invisible(quantiseqr::run_quantiseq(expression_data = tpm,
                                      signature_matrix = "TIL10",
                                      is_tumordata = TRUE))
})
add_bench("quanTIseq", t, paste(ncol(tpm), "amostras"))

t <- system.time({
  invisible(xCell::xCellAnalysis(X_meta, rnaseq = TRUE, parallel.sz = 1))
})
add_bench("xCell", t, paste(ncol(X_meta), "amostras"))

# -----------------------------------------------------------------------------
# Consolidação
# -----------------------------------------------------------------------------
bench <- do.call(rbind, bench)
bench$memoria_max_MB <- round(max(gc()[, 6]), 1)
bench$tempo_total_segundos <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)
fwrite(bench, file.path(dir_aud, "benchmarks_analises_avancadas.csv"))

# DEGs e comparações
resumo_modelos <- fread(file.path(dir_out, "tabelas",
                                   "expressao_diferencial_covariaveis",
                                   "resumo_modelos_expressao_diferencial_covariaveis.csv"))
resumo_combat <- fread(file.path(dir_out, "tabelas", "resumo_testes_ComBat.csv"))
resumo_comp <- fread(file.path(dir_out, "tabelas", "integracao_composicao",
                                "comparacao_base_vs_composicao.csv"))

fwrite(resumo_modelos, file.path(dir_aud, "benchmarks_modelos_expressao_diferencial.csv"))
fwrite(resumo_combat, file.path(dir_aud, "benchmarks_testes_ComBat.csv"))
fwrite(resumo_comp, file.path(dir_aud, "benchmarks_integracao_composicao.csv"))

# Markdown
md <- c(
  "# Benchmarks e consolidação — análises avançadas",
  "",
  paste("Gerado em", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## Tempo e memória",
  "",
  knitr::kable(bench, format = "pipe", digits = 2),
  "",
  "## Modelos de expressão diferencial",
  "",
  knitr::kable(resumo_modelos[, c("modelo", "n_total", "n_covariaveis",
                                  "up", "down", "total_degs", "status")],
               format = "pipe"),
  "",
  "## Testes de ComBat",
  "",
  knitr::kable(resumo_combat[, c("cenario", "status", "n", "n_lotes",
                                 "degs_antes", "degs_depois", "jaccard",
                                 "spearman", "var_lote_pc_antes",
                                 "var_lote_pc_depois", "permanova_antes",
                                 "permanova_depois")],
               format = "pipe", digits = 3),
  "",
  "## Base vs. composição celular",
  "",
  knitr::kable(resumo_comp, format = "pipe", digits = 3),
  "",
  "## Notas",
  "",
  "- Tempos medidos no mesmo computador, sem paralelização do ComBat.",
  "- O tempo do xCell inclui a preparação das assinaturas e o ssGSEA.",
  "- O uso de memória é o máximo observado pelo R durante esta execução.",
  "- Jaccard = interseção/união dos DEGs; Spearman = correlação dos logFC.",
  "- Valores de Jaccard próximos de 0 e Spearman baixo indicam que a correção",
  "  ou o ajuste mudou substancialmente o conjunto de genes.",
  ""
)
writeLines(md, file.path(dir_aud, "benchmarks_analises_avancadas.md"))

cat("Benchmark concluído.\n")
print(bench[, c("etapa", "tempo_segundos", "detalhe")])
