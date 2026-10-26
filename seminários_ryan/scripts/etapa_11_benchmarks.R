# =============================================================================
# etapa_11_benchmarks.R
# -----------------------------------------------------------------------------
# Benchmarks computacionais e de concordância estatística.
#
# Saídas:
#   auditorias/benchmarks_pipeline.md
#   auditorias/benchmarks_pipeline_tabela.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(knitr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_proc <- file.path(ROOT, "dados", "processed")
dir_de   <- file.path(ROOT, "resultados", "tabelas", "expressao_diferencial")
dir_aud  <- file.path(ROOT, "auditorias")
dir.create(dir_aud, recursive = TRUE, showWarnings = FALSE)

FDR_CUT <- 0.05
LFC_CUT <- 1

# -----------------------------------------------------------------------------
# 0. Utilitários
# -----------------------------------------------------------------------------
jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

bench <- list()
add_bench <- function(nome, valor, unidade, observacao = "") {
  bench[[length(bench) + 1]] <<- data.frame(
    benchmark = nome, valor = valor, unidade = unidade,
    observacao = observacao, stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1. Tempo de carga
# -----------------------------------------------------------------------------
t_load <- system.time({
  meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
                sep = "\t", data.table = FALSE, na.strings = c("", "NA"))
  expr_dt <- fread(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"),
                   sep = "\t", data.table = FALSE)
  rownames(expr_dt) <- expr_dt$gene
  expr <- as.matrix(expr_dt[, setdiff(names(expr_dt), "gene")])
  storage.mode(expr) <- "numeric"
})
add_bench("Carga de metadados + matriz", round(t_load[["elapsed"]], 3), "s",
          paste0(nrow(expr), " genes x ", ncol(expr), " amostras"))
add_bench("Memória da matriz", round(as.numeric(object.size(expr)) / 1e6, 1), "MB",
          "matriz numérica em memória")

# -----------------------------------------------------------------------------
# 2. Tempo do UMAP
# -----------------------------------------------------------------------------
vars <- apply(expr, 1, var, na.rm = TRUE)
topv <- names(sort(vars, decreasing = TRUE))[seq_len(min(2000, length(vars)))]
set.seed(2026)
t_umap <- system.time({
  invisible(uwot::umap(t(expr[topv, ]), n_neighbors = 30, min_dist = 0.25,
                       metric = "cosine", scale = TRUE, pca = 50,
                       seed = 2026, n_threads = 1, n_sgd_threads = 1,
                       verbose = FALSE))
})
add_bench("UMAP (2.000 genes, 529 amostras)", round(t_umap[["elapsed"]], 3), "s",
          "n_neighbors=30, min_dist=0.25, cosine, pca=50")

# -----------------------------------------------------------------------------
# 3. Tempo do limma (tumor TCGA vs GTEx normal)
# -----------------------------------------------------------------------------
meta_tg <- meta[meta$group %in% c("TCGA_Tumor", "GTEx_Normal"), ]
meta_tg <- meta_tg[complete.cases(meta_tg[, c("sex", "age")]), ]
meta_tg$group_tg <- factor(ifelse(meta_tg$group == "TCGA_Tumor",
                                  "TCGA_Tumor", "GTEx_Normal"),
                           levels = c("GTEx_Normal", "TCGA_Tumor"))
X <- expr[, meta_tg$sample, drop = FALSE]
design <- model.matrix(~ 0 + group_tg + sex + age, data = meta_tg)
colnames(design) <- sub("group_tg", "group", colnames(design), fixed = TRUE)
t_limma <- system.time({
  fit <- lmFit(X, design)
  cvec <- rep(0, ncol(design)); names(cvec) <- colnames(design)
  cvec["groupGTEx_Normal"] <- -1; cvec["groupTCGA_Tumor"] <- 1
  fit2 <- contrasts.fit(fit, cvec)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  tt_cov <- topTable(fit2, coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt_cov$gene <- rownames(tt_cov)
})
add_bench("limma com covariáveis (tumor vs GTEx)", round(t_limma[["elapsed"]], 3), "s",
          paste0(nrow(X), " genes x ", ncol(X), " amostras"))

# Sem covariáveis (mesmo contraste)
design0 <- model.matrix(~ 0 + group_tg, data = meta_tg)
colnames(design0) <- c("GTEx_Normal", "TCGA_Tumor")
t_limma0 <- system.time({
  fit0 <- lmFit(X, design0)
  fit20 <- contrasts.fit(fit0, c(-1, 1))
  fit20 <- eBayes(fit20, trend = TRUE, robust = TRUE)
  tt_nocov <- topTable(fit20, coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt_nocov$gene <- rownames(tt_nocov)
})
add_bench("limma sem covariáveis (tumor vs GTEx)", round(t_limma0[["elapsed"]], 3), "s",
          "mesmo contraste, sem sexo/idade")

# Trend FALSE (sensibilidade)
t_limma_notrend <- system.time({
  fit_nt <- lmFit(X, design)
  fit2_nt <- contrasts.fit(fit_nt, cvec)
  fit2_nt <- eBayes(fit2_nt, trend = FALSE, robust = TRUE)
  tt_notrend <- topTable(fit2_nt, coef = 1, number = Inf,
                         sort.by = "P", adjust.method = "BH")
  tt_notrend$gene <- rownames(tt_notrend)
})
add_bench("limma com trend = FALSE", round(t_limma_notrend[["elapsed"]], 3), "s",
          "sensibilidade do modelo de variância")

# -----------------------------------------------------------------------------
# 4. Concordância estatística
# -----------------------------------------------------------------------------
sig <- function(d) d$gene[d$adj.P.Val < FDR_CUT & abs(d$logFC) > LFC_CUT]
fdr_only <- function(d) d$gene[d$adj.P.Val < FDR_CUT]

# Com vs sem covariáveis
s_cov <- sig(tt_cov); s_nocov <- sig(tt_nocov)
add_bench("DEGs com covariáveis", length(s_cov), "genes", "FDR<0.05 e |logFC|>1")
add_bench("DEGs sem covariáveis", length(s_nocov), "genes", "FDR<0.05 e |logFC|>1")
add_bench("Jaccard com vs sem covariáveis", round(jaccard(s_cov, s_nocov), 3), "índice",
          "interseção/união dos DEGs")
common <- intersect(tt_cov$gene, tt_nocov$gene)
add_bench("Correlação de logFC (com vs sem covariáveis)",
          round(cor(tt_cov$logFC[match(common, tt_cov$gene)],
                    tt_nocov$logFC[match(common, tt_nocov$gene)],
                    method = "spearman"), 3), "rho", "Spearman")

# trend TRUE vs FALSE
s_trend <- sig(tt_cov); s_notrend <- sig(tt_notrend)
add_bench("DEGs trend = TRUE", length(s_trend), "genes", "")
add_bench("DEGs trend = FALSE", length(s_notrend), "genes", "")
add_bench("Jaccard trend TRUE vs FALSE", round(jaccard(s_trend, s_notrend), 3), "índice", "")

# FDR apenas vs FDR + logFC
s_fdr <- fdr_only(tt_cov)
add_bench("Genes significativos apenas por FDR", length(s_fdr), "genes", "FDR<0.05")
add_bench("Genes significativos por FDR + logFC", length(s_trend), "genes",
          "FDR<0.05 e |logFC|>1")
add_bench("Jaccard FDR vs FDR+logFC", round(jaccard(s_fdr, s_trend), 3), "índice", "")

# Paired vs unpaired
unp <- fread(file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_TCGA_Adjacent_completo.csv"))
par <- fread(file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_Adjacent_PAREADO_completo.csv"))
s_unp <- sig(unp); s_par <- sig(par)
add_bench("DEGs tumor vs adjacente (não pareado)", length(s_unp), "genes", "")
add_bench("DEGs tumor vs adjacente (pareado)", length(s_par), "genes", "")
add_bench("Jaccard pareado vs não pareado", round(jaccard(s_unp, s_par), 3), "índice", "")
common2 <- intersect(unp$gene, par$gene)
add_bench("Correlação de logFC (pareado vs não pareado)",
          round(cor(unp$logFC[match(common2, unp$gene)],
                    par$logFC[match(common2, par$gene)],
                    method = "spearman"), 3), "rho", "Spearman")

# ComBat antes vs depois
antes <- fread(file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_completo.csv"))
depois <- fread(file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_ComBat_completo.csv"))
s_a <- sig(antes); s_d <- sig(depois)
add_bench("DEGs antes do ComBat", length(s_a), "genes", "")
add_bench("DEGs depois do ComBat (sem covariável)", length(s_d), "genes", "")
add_bench("Jaccard ComBat antes vs depois", round(jaccard(s_a, s_d), 3), "índice", "")
common3 <- intersect(antes$gene, depois$gene)
add_bench("Correlação de logFC (ComBat antes vs depois)",
          round(cor(antes$logFC[match(common3, antes$gene)],
                    depois$logFC[match(common3, depois$gene)],
                    method = "spearman"), 3), "rho", "Spearman")

bench_df <- do.call(rbind, bench)
write.csv(bench_df, file.path(dir_aud, "benchmarks_pipeline_tabela.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Relatório markdown
# -----------------------------------------------------------------------------
md <- c(
  "# Benchmarks computacionais e de concordância",
  "",
  paste0("**Data:** ", Sys.Date()),
  "",
  "## 1. Tempo de execução",
  "",
  knitr::kable(bench_df[bench_df$unidade == "s", ], format = "pipe"),
  "",
  "## 2. Concordância estatística",
  "",
  knitr::kable(bench_df[bench_df$unidade != "s", ], format = "pipe"),
  "",
  "## 3. Interpretação",
  "",
  "- O **limma** é rápido o suficiente para 12.401 genes × 478 amostras; o",
  "  custo dominante é a leitura da matriz e a geração das figuras.",
  "- O **UMAP** é viável em poucos segundos com 2.000 genes e 529 amostras,",
  "  usando `n_threads = 1` para maior determinismo.",
  "- A inclusão de **sexo e idade** altera parte dos DEGs; a correlação de",
  "  logFC entre modelos costuma ser alta, mas o conjunto de genes",
  "  significativos não é idêntico. Isso justifica reportar o modelo ajustado",
  "  como principal.",
  "- A comparação **pareado vs não pareado** mostra alta concordância, o que",
  "  reforça a robustez do sinal tumor vs. adjacente.",
  "- A demonstração do **ComBat** reduz drasticamente os DEGs; o Jaccard baixo",
  "  e a correlação de logFC baixa confirmam a sobrecorreção.",
  "- O critério **FDR + |logFC| > 1** é mais conservador que FDR apenas; os",
  "  genes excluídos pelo corte de logFC podem ser consultados nas tabelas",
  "  completas.",
  "",
  "## 4. Arquivos",
  "",
  "- `auditorias/benchmarks_pipeline_tabela.csv`",
  "- `resultados/tabelas/expressao_diferencial/`",
  "",
  "## 5. Sessão R",
  "",
  "```text",
  capture.output(sessionInfo()),
  "```"
)
writeLines(md, file.path(dir_aud, "benchmarks_pipeline.md"))
cat("Benchmarks gerados em:", dir_aud, "\n")
print(bench_df, row.names = FALSE)
