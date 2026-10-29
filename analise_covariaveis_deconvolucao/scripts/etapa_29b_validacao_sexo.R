#!/usr/bin/env Rscript
# Validação externa do dimorfismo sexual no tumor (GSE14520).
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_val <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "validacao_externa")
dir_de  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "expressao_diferencial_covariaveis")
FDR <- 0.05; LFC <- 1
auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label); score <- score[ok]; label <- as.logical(label[ok])
  n1 <- sum(label); n0 <- sum(!label); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label]) - n1*(n1+1)/2)/(n1*n0)
}
Xv_dt <- fread(file.path(dir_val, "GSE14520_expression.tsv"), data.table = FALSE)
Xv <- as.matrix(Xv_dt[, -1, drop = FALSE]); rownames(Xv) <- Xv_dt[[1]]; storage.mode(Xv) <- "numeric"
mv <- fread(file.path(dir_val, "GSE14520_metadata.tsv"), data.table = FALSE)
mv <- mv[mv$tissue == "Tumor" & mv$sex %in% c("M", "F") & !is.na(mv$age), , drop = FALSE]
mv <- mv[mv$gsm %in% colnames(Xv), , drop = FALSE]
Xv <- Xv[, mv$gsm, drop = FALSE]
mv$sex <- factor(mv$sex, levels = c("F", "M"))
design <- model.matrix(~ 0 + sex + age, data = mv)
colnames(design) <- sub("sex", "group", colnames(design))
fit <- lmFit(Xv, design)
cvec <- setNames(rep(0, ncol(design)), colnames(design)); cvec["groupM"] <- 1; cvec["groupF"] <- -1
tt <- topTable(eBayes(contrasts.fit(fit, cvec), trend = TRUE, robust = TRUE),
               coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
tt$gene <- rownames(tt)
tt$regulation <- ifelse(tt$adj.P.Val < FDR & tt$logFC > LFC, "Up",
                 ifelse(tt$adj.P.Val < FDR & tt$logFC < -LFC, "Down", "NS"))
fwrite(tt, file.path(dir_val, "GSE14520_sexo_tumor.csv"))
cat("GSE14520 sexo:", nrow(mv), "tumores;", sum(tt$regulation != "NS"), "DEGs\n")
tc <- fread(file.path(dir_de, "modelo_sexo_tumor_completo.csv"), data.table = FALSE)
com <- intersect(tc$gene, tt$gene)
a <- tc$gene[tc$regulation != "NS"]; b <- tt$gene[tt$regulation != "NS"]
ov <- intersect(a, b)
sig_up <- tc$gene[tc$regulation == "Up"]; sig_dn <- tc$gene[tc$regulation == "Down"]
up_v <- intersect(sig_up, rownames(Xv)); dn_v <- intersect(sig_dn, rownames(Xv))
Z <- t(scale(t(Xv)))
score <- colMeans(Z[up_v, , drop = FALSE], na.rm = TRUE) - colMeans(Z[dn_v, , drop = FALSE], na.rm = TRUE)
res <- data.frame(
  contraste = "Masculino vs feminino no tumor",
  coorte_treinamento = "TCGA-LIHC (n = 342)", coorte_validacao = sprintf("GSE14520 (n = %d)", nrow(mv)),
  genes_comuns = length(com), degs_tcga = length(a), degs_validacao = length(b), degs_comuns = length(ov),
  jaccard = length(ov)/length(union(a, b)),
  hipergeometrico_p = phyper(length(ov)-1, length(intersect(a, com)), length(com)-length(intersect(a, com)),
                             length(intersect(b, com)), lower.tail = FALSE),
  spearman_logfc = cor(tc$logFC[match(com, tc$gene)], tt$logFC[match(com, tt$gene)],
                       method = "spearman", use = "complete.obs"),
  direcao_concordante = mean(sign(tc$logFC[match(ov, tc$gene)]) == sign(tt$logFC[match(ov, tt$gene)])),
  auc_assinatura_tcga = auc(score, mv$sex == "M"),
  xist_p = tt$P.Value[match("XIST", tt$gene)],
  rps4y1_p = tt$P.Value[match("RPS4Y1", tt$gene)],
  stringsAsFactors = FALSE)
fwrite(res, file.path(dir_val, "validacao_sexo.csv"))
print(res)
