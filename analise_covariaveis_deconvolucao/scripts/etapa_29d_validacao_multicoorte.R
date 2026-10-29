#!/usr/bin/env Rscript
# =============================================================================
# etapa_29d_validacao_multicoorte.R
# Validação em coortes independentes públicas: GSE14520, GSE36376, GSE76427
# e GSE62232. Não é uma validação externa prospectiva; são estudos públicos
# independentes, com plataformas e populações diferentes.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(limma); library(dplyr)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_val <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "validacao_externa")
dir_de  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "expressao_diferencial_covariaveis")
FDR <- 0.05; LFC <- 1

auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]; label <- as.logical(label[ok])
  n1 <- sum(label); n0 <- sum(!label)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
carregar <- function(gse) {
  x <- fread(file.path(dir_val, paste0(gse, "_expression.tsv")), data.table = FALSE)
  X <- as.matrix(x[, -1, drop = FALSE]); rownames(X) <- x[[1]]; storage.mode(X) <- "numeric"
  m <- fread(file.path(dir_val, paste0(gse, "_metadata.tsv")), data.table = FALSE)
  m <- m[m$gsm %in% colnames(X), , drop = FALSE]
  m$sex <- toupper(as.character(m$sex))
  m$sex[m$sex %in% c("1", "MALE")] <- "M"
  m$sex[m$sex %in% c("2", "FEMALE")] <- "F"
  m$age <- suppressWarnings(as.numeric(m$age))
  if (!"grade_group" %in% names(m)) m$grade_group <- NA
  if (!"stage_group" %in% names(m)) m$stage_group <- NA
  m$grade_group <- factor(m$grade_group, levels = c("G1-G2", "G3-G4"))
  m$stage_group <- factor(m$stage_group, levels = c("I-II", "III-IV"))
  X <- X[, m$gsm, drop = FALSE]
  list(X = X, m = m)
}

rodar <- function(X, d, grupo, niveis, covs = character(), block_col = NULL) {
  d <- d[d[[grupo]] %in% niveis & !is.na(d[[grupo]]), , drop = FALSE]
  d[[grupo]] <- factor(as.character(d[[grupo]]), levels = niveis)
  keep <- complete.cases(d[, c(grupo, covs), drop = FALSE])
  d <- d[keep, , drop = FALSE]
  if (nrow(d) < 12) return(NULL)
  for (v in names(d)) if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  if (length(covs)) covs <- covs[vapply(d[covs], function(x) length(unique(na.omit(x))) > 1, logical(1))]
  if (length(unique(na.omit(d[[grupo]]))) < 2) return(NULL)
  fml <- as.formula(paste("~ 0 +", grupo, if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""))
  design <- model.matrix(fml, data = d)
  colnames(design) <- sub(grupo, "group", colnames(design), fixed = TRUE)
  if (qr(design)$rank < ncol(design)) return(NULL)
  blk <- if (!is.null(block_col)) d[[block_col]] else NULL
  if (!is.null(blk) && length(unique(blk[keep])) > 1) {
    dup <- duplicateCorrelation(X[, d$gsm, drop = FALSE], design, block = blk[keep])
    fit <- lmFit(X[, d$gsm, drop = FALSE], design, block = blk[keep], correlation = dup$consensus)
  } else {
    fit <- lmFit(X[, d$gsm, drop = FALSE], design)
  }
  cvec <- setNames(rep(0, ncol(design)), colnames(design))
  cvec[paste0("group", niveis[1])] <- -1
  cvec[paste0("group", niveis[2])] <- 1
  if (!all(names(cvec) %in% colnames(design))) return(NULL)
  tt <- topTable(eBayes(contrasts.fit(fit, cvec), trend = TRUE, robust = TRUE),
                 coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt$regulation <- ifelse(tt$adj.P.Val < FDR & tt$logFC > LFC, "Up",
                   ifelse(tt$adj.P.Val < FDR & tt$logFC < -LFC, "Down", "NS"))
  tt
}

comparar <- function(tt, ref_file, cohort, contraste, n1, n2, X, d, label) {
  ref <- fread(file.path(dir_de, ref_file), data.table = FALSE)
  a <- ref$gene[ref$regulation != "NS"]
  com <- intersect(ref$gene, tt$gene)
  b <- tt$gene[tt$regulation != "NS"]
  ov <- intersect(a, b)
  up <- ref$gene[ref$regulation == "Up"]; dn <- ref$gene[ref$regulation == "Down"]
  up <- intersect(up, tt$gene); dn <- intersect(dn, tt$gene)
  Z <- t(scale(t(X[, d$gsm, drop = FALSE])))
  score <- if (length(up) && length(dn)) colMeans(Z[up, , drop = FALSE], na.rm = TRUE) -
             colMeans(Z[dn, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(d))
  data.frame(
    contraste = contraste, coorte = cohort, n_grupo1 = n1, n_grupo2 = n2,
    genes_comuns = length(com), degs_tcga = length(a), degs_coorte = length(b),
    degs_comuns = length(ov),
    jaccard = if (length(union(a, b))) length(ov) / length(union(a, b)) else NA,
    hipergeometrico_p = if (length(com)) phyper(length(ov) - 1, length(intersect(a, com)),
                                               length(com) - length(intersect(a, com)),
                                               length(intersect(b, com)), lower.tail = FALSE) else NA,
    spearman_logfc = if (length(com) > 3) cor(ref$logFC[match(com, ref$gene)],
                                              tt$logFC[match(com, tt$gene)],
                                              method = "spearman", use = "complete.obs") else NA,
    direcao_concordante = if (length(ov)) mean(sign(ref$logFC[match(ov, ref$gene)]) ==
                                                sign(tt$logFC[match(ov, tt$gene)])) else NA,
    assinatura_genes = length(up) + length(dn),
    auc_tcga = auc(score, label),
    stringsAsFactors = FALSE)
}

res <- list()
add <- function(x) if (!is.null(x)) res[[length(res) + 1]] <<- x

# GSE14520: tumor vs adjacente, estágio, sexo
g <- carregar("GSE14520")
X <- g$X; m <- g$m
tt <- rodar(X, m, "tissue", c("Non-Tumor", "Tumor"), character(), block_col = "pair")
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE14520_tumor_adjacente.csv"))
  add(comparar(tt, "modelo_tumor_vs_adjacente_base_completo.csv", "GSE14520",
               "Tumor vs adjacente", sum(m$tissue == "Non-Tumor"), sum(m$tissue == "Tumor"), X, m, m$tissue == "Tumor"))
}
d <- m[m$tissue == "Tumor", , drop = FALSE]
tt <- rodar(X, d, "stage_group", c("I-II", "III-IV"), c("sex", "age"))
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE14520_estagio.csv"))
  add(comparar(tt, "modelo_estagio_base_completo.csv", "GSE14520", "Estágio III-IV vs I-II",
               sum(d$stage_group == "I-II", na.rm=TRUE), sum(d$stage_group == "III-IV", na.rm=TRUE), X, d, d$stage_group == "III-IV"))
}
tt <- rodar(X, d, "sex", c("F", "M"), c("age"))
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE14520_sexo.csv"))
  add(comparar(tt, "modelo_sexo_tumor_completo.csv", "GSE14520", "Masculino vs feminino",
               sum(d$sex == "F", na.rm=TRUE), sum(d$sex == "M", na.rm=TRUE), X, d, d$sex == "M"))
}

# GSE36376: tumor vs adjacente, estágio
g <- carregar("GSE36376"); X <- g$X; m <- g$m
tt <- rodar(X, m, "tissue", c("Non-Tumor", "Tumor"), character())
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE36376_tumor_adjacente.csv"))
  add(comparar(tt, "modelo_tumor_vs_adjacente_base_completo.csv", "GSE36376",
               "Tumor vs adjacente", sum(m$tissue == "Non-Tumor"), sum(m$tissue == "Tumor"), X, m, m$tissue == "Tumor"))
}
d <- m[m$tissue == "Tumor" & !is.na(m$stage_group), , drop = FALSE]
tt <- rodar(X, d, "stage_group", c("I-II", "III-IV"), character())
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE36376_estagio.csv"))
  add(comparar(tt, "modelo_estagio_base_completo.csv", "GSE36376", "Estágio III-IV vs I-II",
               sum(d$stage_group == "I-II", na.rm=TRUE), sum(d$stage_group == "III-IV", na.rm=TRUE), X, d, d$stage_group == "III-IV"))
}

# GSE76427: tumor vs adjacente (pareado), estágio, sexo
g <- carregar("GSE76427"); X <- g$X; m <- g$m
tt <- rodar(X, m, "tissue", c("Non-Tumor", "Tumor"), character(), block_col = "pair")
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE76427_tumor_adjacente.csv"))
  add(comparar(tt, "modelo_tumor_vs_adjacente_base_completo.csv", "GSE76427",
               "Tumor vs adjacente", sum(m$tissue == "Non-Tumor"), sum(m$tissue == "Tumor"), X, m, m$tissue == "Tumor"))
}
d <- m[m$tissue == "Tumor" & !is.na(m$stage_group), , drop = FALSE]
covs <- if (length(unique(na.omit(d$sex))) > 1) c("sex", "age") else "age"
tt <- rodar(X, d, "stage_group", c("I-II", "III-IV"), covs)
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE76427_estagio.csv"))
  add(comparar(tt, "modelo_estagio_base_completo.csv", "GSE76427", "Estágio III-IV vs I-II",
               sum(d$stage_group == "I-II", na.rm=TRUE), sum(d$stage_group == "III-IV", na.rm=TRUE), X, d, d$stage_group == "III-IV"))
}
tt <- rodar(X, d, "sex", c("F", "M"), "age")
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE76427_sexo.csv"))
  add(comparar(tt, "modelo_sexo_tumor_completo.csv", "GSE76427", "Masculino vs feminino",
               sum(d$sex == "F", na.rm=TRUE), sum(d$sex == "M", na.rm=TRUE), X, d, d$sex == "M"))
}

# GSE62232: grau (Edmondson) e sexo
g <- carregar("GSE62232"); X <- g$X; m <- g$m
d <- m[m$tissue == "Tumor" & !is.na(m$grade_group), , drop = FALSE]
tt <- rodar(X, d, "grade_group", c("G1-G2", "G3-G4"), c("sex", "age"))
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE62232_grau.csv"))
  add(comparar(tt, "modelo_grau_base_completo.csv", "GSE62232", "Grau G3-G4 vs G1-G2",
               sum(d$grade_group == "G1-G2", na.rm=TRUE), sum(d$grade_group == "G3-G4", na.rm=TRUE), X, d, d$grade_group == "G3-G4"))
}
tt <- rodar(X, d, "sex", c("F", "M"), "age")
if (!is.null(tt)) {
  fwrite(tt, file.path(dir_val, "multicoorte_GSE62232_sexo.csv"))
  add(comparar(tt, "modelo_sexo_tumor_completo.csv", "GSE62232", "Masculino vs feminino",
               sum(d$sex == "F", na.rm=TRUE), sum(d$sex == "M", na.rm=TRUE), X, d, d$sex == "M"))
}

out <- do.call(rbind, res)
fwrite(out, file.path(dir_val, "validacao_multicoorte.csv"))
print(out)

# Resumo por contraste
meta <- out %>%
  dplyr::group_by(contraste) %>%
  dplyr::summarise(coortes = dplyr::n(), degs_total = sum(degs_coorte, na.rm = TRUE),
                   comuns_total = sum(degs_comuns, na.rm = TRUE),
                   jaccard_mediano = median(jaccard, na.rm = TRUE),
                   spearman_mediano = median(spearman_logfc, na.rm = TRUE),
                   direcao_mediana = median(direcao_concordante, na.rm = TRUE),
                   .groups = "drop")
fwrite(meta, file.path(dir_val, "validacao_multicoorte_resumo.csv"))
print(meta)
