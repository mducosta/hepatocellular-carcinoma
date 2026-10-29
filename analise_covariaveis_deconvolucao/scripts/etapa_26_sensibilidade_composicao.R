#!/usr/bin/env Rscript
# =============================================================================
# etapa_26_sensibilidade_composicao.R
# Sensibilidade do ajuste por composição celular.
#
# O xCell Hepatocytes é estimado da mesma matriz de expressão usada como
# desfecho. Para testar se o desaparecimento dos DEGs de estágio e grau é
# sobreajuste, os modelos são repetidos em universos reduzidos:
#   (a) exclusão de marcadores hepáticos canônicos;
#   (b) exclusão de genes com |Spearman| >= 0,60 com xCell Hepatocytes;
#   (c) interseção de (a) e (b).
# Também é ajustado um escore de marcadores hepáticos curado, independente do
# xCell, como covariável alternativa.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(limma); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))

dir_proc <- file.path(ROOT, "dados", "processed")
dir_dec  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "deconvolucao")
dir_out  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "sensibilidades")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
FDR_CUT <- 0.05; LFC_CUT <- 1

X <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))
uni <- fread(file.path(dir_dec, "deconvolucao_unificada_amostras.csv"), sep = ",",
             data.table = FALSE, na.strings = c("", "NA"))
uni <- uni[uni$sample %in% colnames(X), , drop = FALSE]
for (v in c("group", "sex", "stage_group", "grade_group", "tss_collapsed"))
  if (v %in% names(uni)) uni[[v]] <- factor(uni[[v]])
uni$age <- as.numeric(uni$age)

# Marcadores hepáticos canônicos (curados, sem usar o xCell)
marcadores <- c("ALB","APOA1","APOA2","APOB","APOC1","APOC2","APOC3","APOE",
  "TTR","TF","HP","HPX","FGA","FGB","FGG","AHSG","ITIH1","ITIH2","ITIH3","ITIH4",
  "GC","RBP4","TFR2","SERPINA1","SERPINC1","SERPIND1","CYP1A2","CYP2A6","CYP2A7",
  "CYP2C8","CYP2C9","CYP2D6","CYP2E1","CYP3A4","CYP3A5","CYP7A1","CYP8B1",
  "SLC10A1","SLCO1B1","SLCO1B3","GNMT","SDS","ACOT12","GLS2","HNF4A","HNF1A",
  "FOXA1","FOXA2","FOXA3","ONECUT1","CEBPA","CEBPB","TAT","ARG1","ASS1","ASL",
  "OTC","CPS1","G6PC","PCK1","PCK2","FBP1","ALDOB","KHK","SORD","ADH1B","ADH4",
  "ALDH2","GSTA1","GSTA2","UGT2B7","ABCB11","ABCG5","ABCG8","MARC1","MARC2")
marcadores <- intersect(marcadores, rownames(X))

# Correlação de cada gene com xCell Hepatocytes no tumor
d_tum <- uni[uni$group == "TCGA_Tumor", , drop = FALSE]
Xt <- X[, d_tum$sample, drop = FALSE]
hep <- as.numeric(scale(d_tum[["xCell_Hepatocytes"]]))
rho <- suppressWarnings(cor(t(Xt), hep, method = "spearman", use = "pairwise.complete.obs")[, 1])
rho <- rho[is.finite(rho)]
ordem <- order(-abs(rho))
removidos_rho <- names(rho)[abs(rho) >= 0.60]
universo_rho <- setdiff(rownames(X), removidos_rho)
universo_marc <- setdiff(rownames(X), marcadores)
universo_int <- intersect(universo_rho, universo_marc)

fwrite(data.frame(gene = names(rho)[ordem], rho_hepatocytes = rho[ordem],
                  excluido_rho_060 = names(rho)[ordem] %in% removidos_rho),
       file.path(dir_out, "correlacao_genes_xCell_Hepatocytes.csv"))
fwrite(data.frame(marcador = marcadores), file.path(dir_out, "marcadores_hepaticos_usados.csv"))

# Escore curado de marcadores hepáticos (independente do xCell)
Zm <- t(scale(t(X[marcadores, , drop = FALSE])))
escore <- colMeans(Zm, na.rm = TRUE)
uni$z_escore_hepatico <- as.numeric(scale(escore[uni$sample]))

comp <- c("xCell_Hepatocytes", "xCell_ImmuneScore", "xCell_StromaScore")
for (v in comp) uni[[paste0("z_", v)]] <- as.numeric(scale(uni[[v]]))
comp_z <- paste0("z_", comp)

rodar <- function(Xsub, data, grupo, niveis, covs) {
  rodar_limma(Xsub, data, grupo, niveis, covs, fdr_cut = FDR_CUT, lfc_cut = LFC_CUT)
}
res <- list(); det <- list()
for (contraste in c("estagio", "grau")) {
  if (contraste == "estagio") {
    d <- uni[uni$group == "TCGA_Tumor", , drop = FALSE]
    grupo <- "stage_group"; niveis <- c("I-II", "III-IV")
  } else {
    d <- uni[uni$group == "TCGA_Tumor", , drop = FALSE]
    grupo <- "grade_group"; niveis <- c("G1-G2", "G3-G4")
  }
  base <- rodar(X, d, grupo, niveis, c("sex", "age"))
  cenarios <- list(
    list(nome = "base", Xs = X, covs = c("sex", "age"), universo = "completo"),
    list(nome = "xCell", Xs = X, covs = c("sex", "age", comp_z), universo = "completo"),
    list(nome = "xCell_sem_marcadores", Xs = X[universo_marc, , drop = FALSE],
         covs = c("sex", "age", comp_z), universo = "sem_marcadores_canonicos"),
    list(nome = "xCell_sem_correlatos", Xs = X[universo_rho, , drop = FALSE],
         covs = c("sex", "age", comp_z), universo = "sem_|rho|>=0.60"),
    list(nome = "xCell_universo_restrito", Xs = X[universo_int, , drop = FALSE],
         covs = c("sex", "age", comp_z), universo = "sem_marcadores_e_correlatos"),
    list(nome = "escore_curado", Xs = X, covs = c("sex", "age", "z_escore_hepatico"),
         universo = "completo"),
    list(nome = "escore_curado_sem_marcadores", Xs = X[universo_marc, , drop = FALSE],
         covs = c("sex", "age", "z_escore_hepatico"), universo = "sem_marcadores_canonicos")
  )
  for (cs in cenarios) {
    tt <- rodar(cs$Xs, d, grupo, niveis, cs$covs)
    if (is.null(tt)) { message("falhou: ", contraste, " ", cs$nome); next }
    met <- metricas_concordancia(base, tt)
    res[[paste(contraste, cs$nome, sep = "__")]] <- data.frame(
      contraste = contraste, cenario = cs$nome, universo = cs$universo,
      n_genes = nrow(cs$Xs), n_amostras = sum(complete.cases(d[, c(grupo, cs$covs)])),
      degs = met$degs_depois, up = met$genes_depois_up, down = met$genes_depois_down,
      jaccard_vs_base = met$jaccard, spearman_vs_base = met$spearman_logfc,
      stringsAsFactors = FALSE)
    top <- tt[order(tt$P.Value), c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","regulation")]
    det[[paste(contraste, cs$nome, sep = "__")]] <- head(top, 500)
    sig <- tt[tt$regulation != "NS", c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","regulation")]
    fwrite(sig, file.path(dir_out, paste0("degs_", contraste, "__", cs$nome, ".csv")))
  }
}
resumo <- do.call(rbind, res)
fwrite(resumo, file.path(dir_out, "sensibilidade_composicao_resumo.csv"))
for (nm in names(det)) fwrite(det[[nm]], file.path(dir_out, paste0("top500_", nm, ".csv")))
cat("\nResumo:\n"); print(resumo, row.names = FALSE)
cat("\nGenes removidos por rho>=0.60:", length(removidos_rho), "de", length(rho), "\n")
cat("Marcadores canônicos presentes:", length(marcadores), "\n")
