# =============================================================================
# 02_expressao_diferencial.R
# -----------------------------------------------------------------------------
# Expressão diferencial com limma.
#
# Comparações principais:
#   A) TCGA vs TCGA   : tumor primário vs tecido adjacente (e análise pareada)
#   B) GTEx vs GTEx   : masculino vs feminino; idade >=60 vs <60
#   C) TCGA vs GTEx   : tumor primário (TCGA) vs fígado normal (GTEx)
#                       tecido adjacente (TCGA) vs fígado normal (GTEx)
#   D) Subgrupos TCGA : estágio I-II vs III-IV; grau G1-G2 vs G3-G4;
#                       sexo; idade >=65 vs <65
#   E) Via KEGG hsa05417 (Lipídios e Aterosclerose): tumor vs normal GTEx
#   F) Sensibilidade com ComBat (correção de lote TCGA x GTEx)
#
# Modelo limma: log2(TPM + 0,001) ~ 0 + grupo + covariáveis
# Significância: FDR (BH) < 0,05 e |logFC| > 1.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sva)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_proc <- file.path(ROOT, "dados", "processed")
dir_tab  <- file.path(ROOT, "resultados", "tabelas", "DE")
dir_fig  <- file.path(ROOT, "resultados", "figuras")
dir_qc   <- file.path(ROOT, "resultados", "qc")
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

FDR_CUT <- 0.05
LFC_CUT <- 1

# -----------------------------------------------------------------------------
# 0. Dados
# -----------------------------------------------------------------------------
meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
              sep = "\t", data.table = FALSE)
meta$condition <- factor(meta$condition, levels = c("Tumor", "Adjacent", "Normal"))
meta$group <- factor(meta$group, levels = c("TCGA_Tumor", "TCGA_Adjacent", "GTEx_Normal"))
meta$sex <- factor(meta$sex, levels = c("female", "male"))
meta$stage <- factor(meta$stage, levels = c("I", "II", "III", "IV"))
meta$grade <- factor(meta$grade, levels = c("G1", "G2", "G3", "G4"))

expr_file <- file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz")
if (!file.exists(expr_file)) {
  expr_file <- file.path(dir_proc, "expressao_LA_pathway.tsv.gz")
  warning("Matriz completa ausente; usando painel LA.")
}
expr_dt <- fread(expr_file, sep = "\t", data.table = FALSE)
rownames(expr_dt) <- expr_dt$gene
expr <- as.matrix(expr_dt[, setdiff(names(expr_dt), "gene")])
storage.mode(expr) <- "numeric"
rm(expr_dt)

# -----------------------------------------------------------------------------
# 1. Função principal de DE
# -----------------------------------------------------------------------------
run_limma <- function(expr, meta, comparison, group_col, levels, contrast,
                      covariates = NULL, paired = FALSE, subset_expr = TRUE) {
  # levels: vetor com os níveis do grupo na ordem (referência, comparação)
  # contrast: expressão do contraste, ex. "groupTCGA_Tumor - groupGTEx_Normal"
  meta_sub <- meta[meta[[group_col]] %in% levels, , drop = FALSE]
  meta_sub[[group_col]] <- factor(as.character(meta_sub[[group_col]]), levels = levels)

  # Remover amostras com covariáveis ausentes
  if (!is.null(covariates) && length(covariates) > 0) {
    keep <- complete.cases(meta_sub[, covariates, drop = FALSE])
    meta_sub <- meta_sub[keep, , drop = FALSE]
  }

  # Amostras presentes na matriz
  meta_sub <- meta_sub[meta_sub$sample %in% colnames(expr), , drop = FALSE]
  if (nrow(meta_sub) < 6) {
    warning("Poucas amostras em ", comparison, "; ignorando.")
    return(NULL)
  }

  X <- expr[, meta_sub$sample, drop = FALSE]

  if (paired) {
    # Análise pareada: paciente como fator de bloqueio
    meta_sub$patient <- sub("^(TCGA-[^-]+-[^-]+).*$", "\\1", meta_sub$sample)
    if (length(unique(meta_sub$patient)) < nrow(meta_sub)) {
      design <- model.matrix(~ 0 + meta_sub[[group_col]] + factor(meta_sub$patient))
      colnames(design)[seq_along(levels)] <- paste0("group", levels)
      # O contraste é montado apenas nos coeficientes de grupo
      fit <- lmFit(X, design)
      cvec <- rep(0, ncol(design))
      names(cvec) <- colnames(design)
      cvec[grep("^group", names(cvec))[1]] <- -1
      cvec[grep("^group", names(cvec))[2]] <-  1
      fit2 <- contrasts.fit(fit, cvec)
    } else {
      design <- model.matrix(~ 0 + meta_sub[[group_col]])
      colnames(design) <- paste0("group", levels)
      fit <- lmFit(X, design)
      cvec <- rep(0, ncol(design)); names(cvec) <- colnames(design)
      cvec[paste0("group", levels[1])] <- -1
      cvec[paste0("group", levels[2])] <-  1
      fit2 <- contrasts.fit(fit, cvec)
    }
  } else {
    if (!is.null(covariates) && length(covariates) > 0) {
      fml <- as.formula(paste("~ 0 +", group_col, "+",
                              paste(covariates, collapse = " + ")))
      design <- model.matrix(fml, data = meta_sub)
    } else {
      design <- model.matrix(as.formula(paste("~ 0 +", group_col)), data = meta_sub)
    }
    colnames(design) <- sub(group_col, "group", colnames(design), fixed = TRUE)
    fit <- lmFit(X, design)
    cvec <- rep(0, ncol(design)); names(cvec) <- colnames(design)
    cvec[paste0("group", levels[1])] <- -1
    cvec[paste0("group", levels[2])] <-  1
    fit2 <- contrasts.fit(fit, cvec)
  }

  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "P",
                 adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt$comparison <- comparison
  tt$direction <- ifelse(tt$logFC > 0, paste0(levels[2], " > ", levels[1]),
                         paste0(levels[1], " > ", levels[2]))
  tt$regulation <- ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC >  LFC_CUT, "Up",
                   ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC < -LFC_CUT, "Down", "NS"))
  tt <- tt[, c("gene", "comparison", "direction", "logFC", "AveExpr",
               "t", "P.Value", "adj.P.Val", "B", "regulation")]

  # Salvar resultados completos e significativos
  fname <- gsub("[^A-Za-z0-9_]+", "_", comparison)
  write.csv(tt, file.path(dir_tab, paste0("DEG_", fname, "_completo.csv")),
            row.names = FALSE)
  sig <- tt[tt$regulation != "NS", ]
  write.csv(sig, file.path(dir_tab, paste0("DEG_", fname, "_significativos.csv")),
            row.names = FALSE)

  # Top 20 up e down (por FDR e |logFC|)
  top_up <- tt[tt$regulation == "Up", ]
  top_up <- top_up[order(top_up$adj.P.Val, -abs(top_up$logFC)), ][seq_len(min(20, nrow(top_up))), ]
  top_down <- tt[tt$regulation == "Down", ]
  top_down <- top_down[order(top_down$adj.P.Val, -abs(top_down$logFC)), ][seq_len(min(20, nrow(top_down))), ]
  write.csv(rbind(top_up, top_down),
            file.path(dir_tab, paste0("TOP20_", fname, ".csv")), row.names = FALSE)

  n_up <- sum(tt$regulation == "Up")
  n_down <- sum(tt$regulation == "Down")
  cat(sprintf("[DE] %-55s n=%4d | Up=%4d | Down=%4d\n",
              comparison, nrow(tt), n_up, n_down))
  attr(tt, "n_up") <- n_up
  attr(tt, "n_down") <- n_down
  attr(tt, "n_tested") <- nrow(tt)
  tt
}

resultados <- list()

# -----------------------------------------------------------------------------
# 2. Comparações TCGA x TCGA
# -----------------------------------------------------------------------------
meta$group_adj <- ifelse(meta$group == "TCGA_Tumor", "Tumor",
                  ifelse(meta$group == "TCGA_Adjacent", "Adjacent", NA))
meta$group_adj <- factor(meta$group_adj, levels = c("Adjacent", "Tumor"))

res_tcga <- run_limma(
  expr, meta,
  comparison = "TCGA_Tumor_vs_TCGA_Adjacent",
  group_col = "group_adj", levels = c("Adjacent", "Tumor"),
  contrast = "Tumor - Adjacent",
  covariates = c("sex", "age")
)
resultados[["TCGA_Tumor_vs_TCGA_Adjacent"]] <- res_tcga

# Análise pareada (apenas pacientes com tumor e adjacente)
pares <- meta[meta$group %in% c("TCGA_Tumor", "TCGA_Adjacent"), ]
pares$patient <- sub("^(TCGA-[^-]+-[^-]+).*$", "\\1", pares$sample)
tab_pares <- table(pares$patient)
pares_ok <- names(tab_pares)[tab_pares == 2]
if (length(pares_ok) >= 10) {
  meta_pares <- pares[pares$patient %in% pares_ok, ]
  meta_pares$group_adj <- factor(ifelse(meta_pares$group == "TCGA_Tumor", "Tumor", "Adjacent"),
                                 levels = c("Adjacent", "Tumor"))
  # Para a análise pareada, usamos a função com paired=TRUE
  res_pares <- run_limma(
    expr, meta_pares,
    comparison = "TCGA_Tumor_vs_Adjacent_PAREADO",
    group_col = "group_adj", levels = c("Adjacent", "Tumor"),
    contrast = "Tumor - Adjacent",
    covariates = NULL, paired = TRUE
  )
  resultados[["TCGA_Tumor_vs_Adjacent_PAREADO"]] <- res_pares
} else {
  cat("[DE] Poucos pares tumor/adjacente; análise pareada não executada.\n")
}

# -----------------------------------------------------------------------------
# 3. Comparações GTEx x GTEx
# -----------------------------------------------------------------------------
gtex <- meta[meta$group == "GTEx_Normal", ]
gtex$sex <- factor(gtex$sex, levels = c("female", "male"))
gtex$age60 <- factor(ifelse(gtex$age >= 60, ">=60", "<60"), levels = c("<60", ">=60"))

res_gtex_sex <- run_limma(
  expr, gtex,
  comparison = "GTEx_Male_vs_Female",
  group_col = "sex", levels = c("female", "male"),
  contrast = "Male - Female",
  covariates = c("age")
)
resultados[["GTEx_Male_vs_Female"]] <- res_gtex_sex

res_gtex_age <- run_limma(
  expr, gtex,
  comparison = "GTEx_Age60plus_vs_lt60",
  group_col = "age60", levels = c("<60", ">=60"),
  contrast = ">=60 - <60",
  covariates = c("sex")
)
resultados[["GTEx_Age60plus_vs_lt60"]] <- res_gtex_age

# -----------------------------------------------------------------------------
# 4. Comparações TCGA x GTEx
# -----------------------------------------------------------------------------
meta_tg <- meta[meta$group %in% c("TCGA_Tumor", "GTEx_Normal"), ]
meta_tg$group_tg <- factor(ifelse(meta_tg$group == "TCGA_Tumor", "TCGA_Tumor", "GTEx_Normal"),
                           levels = c("GTEx_Normal", "TCGA_Tumor"))

res_tg <- run_limma(
  expr, meta_tg,
  comparison = "TCGA_Tumor_vs_GTEx_Normal",
  group_col = "group_tg", levels = c("GTEx_Normal", "TCGA_Tumor"),
  contrast = "TCGA tumor - GTEx normal",
  covariates = c("sex", "age")
)
resultados[["TCGA_Tumor_vs_GTEx_Normal"]] <- res_tg

# TCGA adjacente vs GTEx normal (avalia se o "normal" TCGA é comparável ao GTEx)
meta_ag <- meta[meta$group %in% c("TCGA_Adjacent", "GTEx_Normal"), ]
meta_ag$group_ag <- factor(ifelse(meta_ag$group == "TCGA_Adjacent", "TCGA_Adjacent", "GTEx_Normal"),
                           levels = c("GTEx_Normal", "TCGA_Adjacent"))

res_ag <- run_limma(
  expr, meta_ag,
  comparison = "TCGA_Adjacent_vs_GTEx_Normal",
  group_col = "group_ag", levels = c("GTEx_Normal", "TCGA_Adjacent"),
  contrast = "TCGA adjacent - GTEx normal",
  covariates = c("sex", "age")
)
resultados[["TCGA_Adjacent_vs_GTEx_Normal"]] <- res_ag

# -----------------------------------------------------------------------------
# 5. Subgrupos clínicos TCGA (apenas tumores)
# -----------------------------------------------------------------------------
tcga <- meta[meta$group == "TCGA_Tumor", ]
tcga$stage2 <- factor(ifelse(tcga$stage %in% c("I", "II"), "I-II",
                      ifelse(tcga$stage %in% c("III", "IV"), "III-IV", NA)),
                      levels = c("I-II", "III-IV"))
tcga$grade2 <- factor(ifelse(tcga$grade %in% c("G1", "G2"), "G1-G2",
                      ifelse(tcga$grade %in% c("G3", "G4"), "G3-G4", NA)),
                      levels = c("G1-G2", "G3-G4"))
tcga$age65 <- factor(ifelse(tcga$age >= 65, ">=65", "<65"), levels = c("<65", ">=65"))
tcga$sex <- factor(tcga$sex, levels = c("female", "male"))

res_stage <- run_limma(
  expr, tcga,
  comparison = "TCGA_Stage_III_IV_vs_I_II",
  group_col = "stage2", levels = c("I-II", "III-IV"),
  contrast = "III-IV - I-II",
  covariates = c("sex", "age")
)
resultados[["TCGA_Stage_III_IV_vs_I_II"]] <- res_stage

res_grade <- run_limma(
  expr, tcga,
  comparison = "TCGA_Grade_G3_G4_vs_G1_G2",
  group_col = "grade2", levels = c("G1-G2", "G3-G4"),
  contrast = "G3-G4 - G1-G2",
  covariates = c("sex", "age")
)
resultados[["TCGA_Grade_G3_G4_vs_G1_G2"]] <- res_grade

# Sexo em TCGA (ajustando por idade)
res_sex_tcga <- run_limma(
  expr, tcga,
  comparison = "TCGA_Male_vs_Female",
  group_col = "sex", levels = c("female", "male"),
  contrast = "Male - Female",
  covariates = c("age")
)
resultados[["TCGA_Male_vs_Female"]] <- res_sex_tcga

res_age_tcga <- run_limma(
  expr, tcga,
  comparison = "TCGA_Age65plus_vs_lt65",
  group_col = "age65", levels = c("<65", ">=65"),
  contrast = ">=65 - <65",
  covariates = c("sex")
)
resultados[["TCGA_Age65plus_vs_lt65"]] <- res_age_tcga

# -----------------------------------------------------------------------------
# 6. Via LA (KEGG hsa05417) — análise específica
# -----------------------------------------------------------------------------
la_file <- file.path(dir_proc, "expressao_LA_pathway.tsv.gz")
if (file.exists(la_file)) {
  la_dt <- fread(la_file, sep = "\t", data.table = FALSE)
  rownames(la_dt) <- la_dt$gene
  la_expr <- as.matrix(la_dt[, setdiff(names(la_dt), "gene")])
  storage.mode(la_expr) <- "numeric"

  res_la <- run_limma(
    la_expr, meta_tg,
    comparison = "LApathway_TCGA_Tumor_vs_GTEx_Normal",
    group_col = "group_tg", levels = c("GTEx_Normal", "TCGA_Tumor"),
    contrast = "TCGA tumor - GTEx normal",
    covariates = c("sex", "age")
  )
  resultados[["LApathway_TCGA_Tumor_vs_GTEx_Normal"]] <- res_la

  res_la_adj <- run_limma(
    la_expr, meta,
    comparison = "LApathway_TCGA_Tumor_vs_TCGA_Adjacent",
    group_col = "group_adj", levels = c("Adjacent", "Tumor"),
    contrast = "Tumor - Adjacent",
    covariates = c("sex", "age")
  )
  resultados[["LApathway_TCGA_Tumor_vs_TCGA_Adjacent"]] <- res_la_adj
}

# -----------------------------------------------------------------------------
# 7. Sensibilidade: ComBat (TCGA x GTEx)
# -----------------------------------------------------------------------------
cat("\n[COMBAT] Corrigindo lote TCGA x GTEx e reavaliando tumor vs normal ...\n")
meta_cb <- meta[meta$group %in% c("TCGA_Tumor", "GTEx_Normal"), ]
meta_cb$group_tg <- factor(ifelse(meta_cb$group == "TCGA_Tumor", "TCGA_Tumor", "GTEx_Normal"),
                           levels = c("GTEx_Normal", "TCGA_Tumor"))
meta_cb <- meta_cb[complete.cases(meta_cb[, c("sex", "age")]), ]
Xcb <- expr[, meta_cb$sample, drop = FALSE]

mod_cb <- model.matrix(~ group_tg, data = meta_cb)
batch_cb <- factor(meta_cb$study, levels = c("TCGA", "GTEX"))

combat_ok <- TRUE
Xcb_cor <- tryCatch(
  sva::ComBat(dat = Xcb, batch = batch_cb, mod = mod_cb, par.prior = TRUE),
  error = function(e) { combat_ok <<- FALSE; message("ComBat falhou: ", e$message); NULL }
)

if (combat_ok && !is.null(Xcb_cor)) {
  res_cb <- run_limma(
    Xcb_cor, meta_cb,
    comparison = "TCGA_Tumor_vs_GTEx_Normal_ComBat",
    group_col = "group_tg", levels = c("GTEx_Normal", "TCGA_Tumor"),
    contrast = "TCGA tumor - GTEx normal",
    covariates = c("sex", "age")
  )
  resultados[["TCGA_Tumor_vs_GTEx_Normal_ComBat"]] <- res_cb

  # UMAP antes/depois do ComBat
  vars <- apply(Xcb, 1, var, na.rm = TRUE)
  topv <- names(sort(vars, decreasing = TRUE))[seq_len(min(2000, length(vars)))]
  set.seed(2026)
  umap_before <- uwot::umap(t(Xcb[topv, ]), n_neighbors = 30, min_dist = 0.25,
                            metric = "cosine", scale = TRUE, pca = 50,
                            seed = 2026, n_threads = 1, n_sgd_threads = 1,
                            verbose = FALSE)
  umap_after <- uwot::umap(t(Xcb_cor[topv, ]), n_neighbors = 30, min_dist = 0.25,
                           metric = "cosine", scale = TRUE, pca = 50,
                           seed = 2026, n_threads = 1, n_sgd_threads = 1,
                           verbose = FALSE)
  umap_cb <- rbind(
    data.frame(sample = colnames(Xcb), UMAP1 = umap_before[, 1],
               UMAP2 = umap_before[, 2], etapa = "Antes do ComBat"),
    data.frame(sample = colnames(Xcb_cor), UMAP1 = umap_after[, 1],
               UMAP2 = umap_after[, 2], etapa = "Depois do ComBat")
  )
  umap_cb <- merge(umap_cb, meta_cb[, c("sample", "group_tg", "study")], by = "sample")
  write.csv(umap_cb, file.path(dir_qc, "UMAP_ComBat_coordenadas.csv"), row.names = FALSE)

  p_cb <- ggplot(umap_cb, aes(UMAP1, UMAP2, colour = group_tg, shape = study)) +
    geom_point(size = 1.7, alpha = 0.75) +
    facet_wrap(~ etapa, scales = "free") +
    scale_colour_manual(values = c("TCGA_Tumor" = "#C0392B",
                                   "GTEx_Normal" = "#2471A3")) +
    labs(title = "Efeito do ComBat na estrutura global (TCGA x GTEx)",
         subtitle = "UMAP com os 2.000 genes mais variáveis, antes e depois da correção de lote",
         x = "UMAP 1", y = "UMAP 2", colour = "Grupo", shape = "Estudo") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  ggsave(file.path(dir_fig, "13_UMAP_ComBat_antes_depois.png"), p_cb,
         width = 10, height = 5, dpi = 300)

  # Comparação dos DEGs antes/depois
  antes <- resultados[["TCGA_Tumor_vs_GTEx_Normal"]]
  depois <- resultados[["TCGA_Tumor_vs_GTEx_Normal_ComBat"]]
  comp_cb <- data.frame(
    metrica = c("DEGs antes (FDR<0.05, |logFC|>1)",
                "DEGs depois (FDR<0.05, |logFC|>1)",
                "Up antes", "Up depois", "Down antes", "Down depois",
                "Genes em comum (DEGs)", "Jaccard dos DEGs"),
    valor = c(
      sum(antes$regulation != "NS"), sum(depois$regulation != "NS"),
      sum(antes$regulation == "Up"), sum(depois$regulation == "Up"),
      sum(antes$regulation == "Down"), sum(depois$regulation == "Down"),
      length(intersect(antes$gene[antes$regulation != "NS"],
                       depois$gene[depois$regulation != "NS"])),
      round(length(intersect(antes$gene[antes$regulation != "NS"],
                             depois$gene[depois$regulation != "NS"])) /
              length(union(antes$gene[antes$regulation != "NS"],
                           depois$gene[depois$regulation != "NS"])), 3)
    )
  )
  write.csv(comp_cb, file.path(dir_tab, "COMPARACAO_ComBat_antes_depois.csv"),
            row.names = FALSE)
} else {
  cat("[COMBAT] Não foi possível executar o ComBat.\n")
}

# -----------------------------------------------------------------------------
# 8. Resumo geral
# -----------------------------------------------------------------------------
resumo <- data.frame(
  comparacao = names(resultados),
  n_testados = sapply(resultados, function(x) attr(x, "n_tested")),
  n_up = sapply(resultados, function(x) attr(x, "n_up")),
  n_down = sapply(resultados, function(x) attr(x, "n_down")),
  n_total = sapply(resultados, function(x) attr(x, "n_up") + attr(x, "n_down")),
  percentual_sig = round(100 * sapply(resultados, function(x)
    (attr(x, "n_up") + attr(x, "n_down")) / attr(x, "n_tested")), 2),
  row.names = NULL
)
resumo <- resumo[order(resumo$n_total, decreasing = TRUE), ]
write.csv(resumo, file.path(dir_tab, "RESUMO_DEGs_TODAS_COMPARACOES.csv"),
          row.names = FALSE)
writeLines(capture.output(print(resumo, row.names = FALSE)),
           file.path(dir_qc, "resumo_expressao_diferencial.txt"))

cat("\nResumo geral:\n")
print(resumo, row.names = FALSE)
cat("\nSessão R:\n")
print(sessionInfo())
