# =============================================================================
# etapa_14_ComBat_testes_multiplos.R
# -----------------------------------------------------------------------------
# Bateria de testes do ComBat (sva) para efeito de lote.
#
# Cenários:
#   A) TCGA tumor vs. GTEx normal, lote = estudo.
#      A1. mod = ~ grupo + sexo + idade  -> confundimento estrutural (erro).
#      A2. mod = ~ sexo + idade          -> remove lote, mas não protege grupo.
#      A3. mod = NULL                    -> sobrecorreção (remove biologia).
#      A4. mean.only = TRUE              -> remove só efeito aditivo.
#      A5. ref.batch = GTEX              -> lote de referência.
#   B) TCGA tumor, lote = centro (TSS), com preservação de grau ou estágio.
#   C) TCGA tumor + adjacente, lote = TSS, preservando tumor/adjacente.
#
# Métricas:
#   - DEGs antes/depois (FDR < 0,05 e |logFC| > 1);
#   - Jaccard e Spearman de logFC;
#   - R² do lote em PCs (top 2000 genes variáveis);
#   - PERMANOVA (adonis2) do lote;
#   - preservação de genes marcadores.
#
# Observação: ComBat clássico é aplicado em log2(TPM + 0,001). Para contagens
# brutas, o correto seria ComBat-seq; aqui não há contagens.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sva)
  library(limma)
  library(vegan)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))

dir_proc <- file.path(ROOT, "dados", "processed")
dir_out  <- file.path(ROOT, "resultados", "analises_avancadas")
dir_tab  <- file.path(dir_out, "tabelas", "ComBat")
dir_pc   <- file.path(dir_out, "dados_derivados", "PCA_ComBat")
dir_qc   <- file.path(dir_out, "qc")
for (d in c(dir_tab, dir_pc, dir_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
FDR_CUT <- 0.05
LFC_CUT <- 1

prep <- preparar_metadados(file.path(dir_proc, "metadata_analise.tsv"), dedup_gtex = TRUE)
meta <- prep$meta
X <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))

# -----------------------------------------------------------------------------
# Métricas
# -----------------------------------------------------------------------------
genes_variaveis <- function(X, n = 2000) {
  v <- matrixStats::rowVars(X)
  names(sort(v, decreasing = TRUE))[seq_len(min(n, length(v)))]
}

r2_lote_pcs <- function(X, genes, batch, npc = 10) {
  pca <- prcomp(t(X[genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
  r2 <- vapply(seq_len(npc), function(i) {
    summary(lm(pca$x[, i] ~ batch))$r.squared
  }, numeric(1))
  list(r2 = r2, weighted = weighted.mean(r2, pca$sdev[seq_len(npc)]^2),
       pca = pca)
}

permanova_lote <- function(X, genes, batch, permutations = 99) {
  dd <- vegan::vegdist(t(X[genes, , drop = FALSE]), method = "euclidean")
  set.seed(2026)
  fit <- vegan::adonis2(dd ~ batch, data = data.frame(batch = batch),
                        permutations = permutations)
  as.numeric(fit$R2[1])
}

avaliar <- function(X, d, group_col, levels, covariates) {
  rodar_limma(X, d, group_col, levels, covariates,
              fdr_cut = FDR_CUT, lfc_cut = LFC_CUT)
}

# -----------------------------------------------------------------------------
# Executor de cenário
# -----------------------------------------------------------------------------
executar_cenario <- function(nome, data, batch_col, formula_mod = NULL,
                             group_col, levels, eval_covariates = character(),
                             mean.only = FALSE, ref.batch = NULL,
                             note = "") {
  message("Cenário: ", nome)

  vars_mod <- if (is.null(formula_mod)) character() else all.vars(formula_mod)
  d0 <- data[complete.cases(data[, unique(c(batch_col, vars_mod)), drop = FALSE]), ,
             drop = FALSE]
  d0[[batch_col]] <- factor(as.character(d0[[batch_col]]))
  d0[[batch_col]] <- droplevels(d0[[batch_col]])

  # Remove lotes com menos de 2 amostras (ComBat não é estável)
  tab_batch <- table(d0[[batch_col]])
  b_ok <- names(tab_batch)[tab_batch >= 2]
  d0 <- d0[d0[[batch_col]] %in% b_ok, , drop = FALSE]
  d0[[batch_col]] <- droplevels(d0[[batch_col]])

  if (nlevels(d0[[batch_col]]) < 2) {
    return(list(summary = data.frame(
      cenario = nome, status = "nao_executado", n = nrow(d0),
      n_lotes = nlevels(d0[[batch_col]]), mod = ifelse(is.null(formula_mod), "NULL",
                                                       deparse(formula_mod)),
      mean_only = mean.only, ref_batch = ifelse(is.null(ref.batch), "nenhum", ref.batch),
      comparacao = paste(levels[2], "vs", levels[1]),
      degs_antes = NA, degs_depois = NA, jaccard = NA, spearman = NA,
      var_lote_pc_antes = NA, var_lote_pc_depois = NA,
      permanova_antes = NA, permanova_depois = NA,
      erro = "menos de dois lotes com n >= 2", observacao = note,
      stringsAsFactors = FALSE), tabela_antes = NULL, tabela_depois = NULL,
      pca = NULL))
  }

  X0 <- X[, d0$sample, drop = FALSE]
  batch0 <- d0[[batch_col]]

  # Tabela antes (sempre que possível)
  d_eval <- d0[complete.cases(d0[, c(group_col, eval_covariates), drop = FALSE]), ,
               drop = FALSE]
  antes <- avaliar(X0[, d_eval$sample, drop = FALSE], d_eval,
                   group_col, levels, eval_covariates)

  # ComBat
  mod <- if (is.null(formula_mod)) NULL else model.matrix(formula_mod, data = d0)
  res_cb <- tryCatch({
    if (is.null(ref.batch)) {
      sva::ComBat(dat = X0, batch = batch0, mod = mod,
                  par.prior = TRUE, mean.only = mean.only)
    } else {
      sva::ComBat(dat = X0, batch = batch0, mod = mod,
                  par.prior = TRUE, mean.only = mean.only, ref.batch = ref.batch)
    }
  }, error = function(e) e)

  status <- "ok"
  erro <- ""
  X_cor <- NULL
  if (inherits(res_cb, "error")) {
    status <- "erro"
    erro <- conditionMessage(res_cb)
  } else {
    X_cor <- res_cb
  }

  depois <- NULL
  if (!is.null(X_cor)) {
    depois <- avaliar(X_cor[, d_eval$sample, drop = FALSE], d_eval,
                      group_col, levels, eval_covariates)
  }
  met <- metricas_concordancia(antes, depois)

  # Métricas de lote na mesma matriz de avaliação
  genes <- genes_variaveis(X0[, d_eval$sample, drop = FALSE], 2000)
  batch_eval <- factor(as.character(d_eval[[batch_col]]))
  pc_antes <- tryCatch(r2_lote_pcs(X0[, d_eval$sample, drop = FALSE], genes, batch_eval),
                       error = function(e) NULL)
  perm_antes <- tryCatch(permanova_lote(X0[, d_eval$sample, drop = FALSE], genes, batch_eval),
                         error = function(e) NA_real_)
  pc_depois <- NULL
  perm_depois <- NA_real_
  if (!is.null(X_cor)) {
    pc_depois <- tryCatch(r2_lote_pcs(X_cor[, d_eval$sample, drop = FALSE], genes, batch_eval),
                          error = function(e) NULL)
    perm_depois <- tryCatch(permanova_lote(X_cor[, d_eval$sample, drop = FALSE], genes, batch_eval),
                            error = function(e) NA_real_)
  }

  # Coordenadas PCA para figuras
  pca_df <- NULL
  if (!is.null(pc_antes)) {
    pca_df <- data.frame(
      sample = d_eval$sample,
      PC1 = pc_antes$pca$x[, 1], PC2 = pc_antes$pca$x[, 2], PC3 = pc_antes$pca$x[, 3],
      batch = as.character(batch_eval),
      group = as.character(d_eval[[group_col]]),
      cenario = nome, momento = "antes", stringsAsFactors = FALSE
    )
    if (!is.null(pc_depois)) {
      pca_df <- rbind(pca_df, data.frame(
        sample = d_eval$sample,
        PC1 = pc_depois$pca$x[, 1], PC2 = pc_depois$pca$x[, 2], PC3 = pc_depois$pca$x[, 3],
        batch = as.character(batch_eval),
        group = as.character(d_eval[[group_col]]),
        cenario = nome, momento = "depois", stringsAsFactors = FALSE
      ))
    }
    fwrite(pca_df, file.path(dir_pc, paste0("PCA_", nome, ".csv")))
  }

  sm <- data.frame(
    cenario = nome, status = status, n = nrow(d0), n_lotes = nlevels(batch0),
    mod = ifelse(is.null(formula_mod), "NULL", paste(deparse(formula_mod), collapse = " ")),
    mean_only = mean.only, ref_batch = ifelse(is.null(ref.batch), "nenhum", ref.batch),
    comparacao = paste(levels[2], "vs", levels[1]),
    degs_antes = met$degs_antes, degs_depois = met$degs_depois,
    jaccard = met$jaccard, spearman = met$spearman_logfc,
    var_lote_pc_antes = ifelse(is.null(pc_antes), NA, pc_antes$weighted),
    var_lote_pc_depois = ifelse(is.null(pc_depois), NA, pc_depois$weighted),
    permanova_antes = perm_antes, permanova_depois = perm_depois,
    erro = erro, observacao = note, stringsAsFactors = FALSE
  )

  list(summary = sm, tabela_antes = antes, tabela_depois = depois, pca = pca_df)
}

# -----------------------------------------------------------------------------
# Cenários
# -----------------------------------------------------------------------------
cenarios <- list(
  list(nome = "tcga_gtex_mod_grupo", data = meta,
       batch_col = "study", formula_mod = ~ group + sex + age,
       group_col = "group", levels = c("GTEx_Normal", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "Tentativa de preservar tumor/normal; confundimento estrutural."),
  list(nome = "tcga_gtex_mod_sexo_idade", data = meta,
       batch_col = "study", formula_mod = ~ sex + age,
       group_col = "group", levels = c("GTEx_Normal", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "Remove lote, mas não protege tumor/normal."),
  list(nome = "tcga_gtex_sem_mod", data = meta,
       batch_col = "study", formula_mod = NULL,
       group_col = "group", levels = c("GTEx_Normal", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "Sem covariável biológica; demonstração de sobrecorreção."),
  list(nome = "tcga_gtex_sem_mod_meanonly", data = meta,
       batch_col = "study", formula_mod = NULL, mean.only = TRUE,
       group_col = "group", levels = c("GTEx_Normal", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "Remove apenas efeito aditivo; sem covariável biológica."),
  list(nome = "tcga_gtex_ref_gtex", data = meta,
       batch_col = "study", formula_mod = NULL, ref.batch = "GTEX",
       group_col = "group", levels = c("GTEx_Normal", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "GTEx como lote de referência; sem covariável biológica."),
  list(nome = "tcga_tss_grau", data = prep$tumor,
       batch_col = "tss_collapsed", formula_mod = ~ grade_group + sex + age,
       group_col = "grade_group", levels = c("G1-G2", "G3-G4"),
       eval_covariates = c("sex", "age"),
       note = "TSS corrigido preservando grau; avaliação sem estágio."),
  list(nome = "tcga_tss_grau_estagio", data = prep$tumor,
       batch_col = "tss_collapsed", formula_mod = ~ grade_group + stage_group + sex + age,
       group_col = "grade_group", levels = c("G1-G2", "G3-G4"),
       eval_covariates = c("sex", "age", "stage_group"),
       note = "TSS corrigido preservando grau e estágio."),
  list(nome = "tcga_tss_grau_estagio_meanonly", data = prep$tumor,
       batch_col = "tss_collapsed", formula_mod = ~ grade_group + stage_group + sex + age,
       mean.only = TRUE,
       group_col = "grade_group", levels = c("G1-G2", "G3-G4"),
       eval_covariates = c("sex", "age", "stage_group"),
       note = "Igual ao anterior, mean.only = TRUE."),
  list(nome = "tcga_tss_grau_refDD", data = prep$tumor,
       batch_col = "tss_collapsed", formula_mod = ~ grade_group + stage_group + sex + age,
       ref.batch = "DD",
       group_col = "grade_group", levels = c("G1-G2", "G3-G4"),
       eval_covariates = c("sex", "age", "stage_group"),
       note = "Centro DD como referência."),
  list(nome = "tcga_tss_estagio", data = prep$tumor,
       batch_col = "tss_collapsed", formula_mod = ~ stage_group + sex + age,
       group_col = "stage_group", levels = c("I-II", "III-IV"),
       eval_covariates = c("sex", "age"),
       note = "TSS corrigido preservando estágio; avaliação sem grau."),
  list(nome = "tcga_tss_estagio_grau", data = prep$tumor,
       batch_col = "tss_collapsed", formula_mod = ~ stage_group + grade_group + sex + age,
       group_col = "stage_group", levels = c("I-II", "III-IV"),
       eval_covariates = c("sex", "age", "grade_group"),
       note = "TSS corrigido preservando estágio e grau."),
  list(nome = "tcga_tss_estagio_grau_meanonly", data = prep$tumor,
       batch_col = "tss_collapsed", formula_mod = ~ stage_group + grade_group + sex + age,
       mean.only = TRUE,
       group_col = "stage_group", levels = c("I-II", "III-IV"),
       eval_covariates = c("sex", "age", "grade_group"),
       note = "Igual ao anterior, mean.only = TRUE."),
  list(nome = "tcga_ta_tss", data = prep$tcga,
       batch_col = "tss_collapsed", formula_mod = ~ group + sex + age,
       group_col = "group", levels = c("TCGA_Adjacent", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "Tumor vs. adjacente com TSS corrigido."),
  list(nome = "tcga_ta_tss_meanonly", data = prep$tcga,
       batch_col = "tss_collapsed", formula_mod = ~ group + sex + age,
       mean.only = TRUE,
       group_col = "group", levels = c("TCGA_Adjacent", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "Tumor vs. adjacente com TSS corrigido, mean.only = TRUE."),
  list(nome = "tcga_ta_tss_mod_sexo_idade", data = prep$tcga,
       batch_col = "tss_collapsed", formula_mod = ~ sex + age,
       group_col = "group", levels = c("TCGA_Adjacent", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "TSS corrigido sem proteger tumor/adjacente; demonstra sobrecorreção."),
  list(nome = "tcga_ta_tss_sem_mod", data = prep$tcga,
       batch_col = "tss_collapsed", formula_mod = NULL,
       group_col = "group", levels = c("TCGA_Adjacent", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "TSS corrigido sem covariável biológica."),
  list(nome = "tcga_ta_tss_sem_mod_meanonly", data = prep$tcga,
       batch_col = "tss_collapsed", formula_mod = NULL, mean.only = TRUE,
       group_col = "group", levels = c("TCGA_Adjacent", "TCGA_Tumor"),
       eval_covariates = c("sex", "age"),
       note = "TSS corrigido sem covariável biológica, mean.only = TRUE.")
)

resultados <- list()
for (cs in cenarios) {
  out <- tryCatch(
    do.call(executar_cenario, cs),
    error = function(e) {
      list(summary = data.frame(
        cenario = cs$nome, status = "erro_funcao", n = NA, n_lotes = NA,
        mod = ifelse(is.null(cs$formula_mod), "NULL", paste(deparse(cs$formula_mod), collapse = " ")),
        mean_only = ifelse(is.null(cs$mean.only), FALSE, cs$mean.only),
        ref_batch = ifelse(is.null(cs$ref.batch), "nenhum", cs$ref.batch),
        comparacao = paste(cs$levels[2], "vs", cs$levels[1]),
        degs_antes = NA, degs_depois = NA, jaccard = NA, spearman = NA,
        var_lote_pc_antes = NA, var_lote_pc_depois = NA,
        permanova_antes = NA, permanova_depois = NA,
        erro = conditionMessage(e), observacao = cs$note, stringsAsFactors = FALSE),
        tabela_antes = NULL, tabela_depois = NULL, pca = NULL)
    }
  )
  resultados[[cs$nome]] <- out
  if (!is.null(out$tabela_antes)) {
    fwrite(out$tabela_antes, file.path(dir_tab, paste0(out$summary$cenario, "_antes.csv")))
  }
  if (!is.null(out$tabela_depois)) {
    fwrite(out$tabela_depois, file.path(dir_tab, paste0(out$summary$cenario, "_depois.csv")))
  }
}

resumo <- do.call(rbind, lapply(resultados, `[[`, "summary"))
rownames(resumo) <- NULL
fwrite(resumo, file.path(dir_out, "tabelas", "resumo_testes_ComBat.csv"))

# Tabela de marcadores antes/depois
marcadores <- c("ALB", "AFP", "GPC3", "BIRC5", "MMP9", "MMP1", "CYP1A2",
                "CLEC4M", "CLEC4G", "XIST", "RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY")
mark_list <- list()
for (nm in names(resultados)) {
  a <- resultados[[nm]]$tabela_antes
  b <- resultados[[nm]]$tabela_depois
  if (is.null(a) || is.null(b)) next
  aa <- a[a$gene %in% marcadores, c("gene", "logFC", "adj.P.Val", "regulation")]
  bb <- b[b$gene %in% marcadores, c("gene", "logFC", "adj.P.Val", "regulation")]
  m <- merge(aa, bb, by = "gene", suffixes = c("_antes", "_depois"))
  m$cenario <- nm
  mark_list[[nm]] <- m
}
if (length(mark_list)) {
  fwrite(do.call(rbind, mark_list), file.path(dir_tab, "marcadores_antes_depois_ComBat.csv"))
}

# -----------------------------------------------------------------------------
# Log
# -----------------------------------------------------------------------------
log_lines <- c(
  "TESTES MÚLTIPLOS DE COMBAT — LOG",
  "================================",
  paste("Data:", Sys.time()),
  paste("R:", R.version.string),
  paste("Genes:", nrow(X), "| Amostras totais:", ncol(X)),
  "",
  paste0(resumo$cenario, ": ", resumo$status,
         " | n=", resumo$n, " | lotes=", resumo$n_lotes,
         " | DEGs antes=", resumo$degs_antes, " depois=", resumo$degs_depois,
         " | Jaccard=", round(resumo$jaccard, 3),
         " | R2 lote antes=", round(resumo$var_lote_pc_antes, 3),
         " depois=", round(resumo$var_lote_pc_depois, 3),
         ifelse(resumo$erro == "", "", paste0(" | erro: ", resumo$erro)))
)
writeLines(log_lines, file.path(dir_qc, "log_testes_multiplos_ComBat.txt"))

cat("Concluído.\n")
print(resumo[, c("cenario", "status", "n", "n_lotes", "degs_antes", "degs_depois",
                 "jaccard", "spearman", "var_lote_pc_antes", "var_lote_pc_depois",
                 "permanova_antes", "permanova_depois")])
