# =============================================================================
# utilitarios_analises_avancadas.R
# -----------------------------------------------------------------------------
# Funções compartilhadas pelas análises avançadas:
#   - preparação e correção dos metadados;
#   - leitura da matriz log2(TPM + 0,001);
#   - conversão para TPM linear;
#   - modelos limma simples;
#   - métricas de concordância entre listas de DEGs.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

preparar_metadados <- function(metadata_path, dedup_gtex = TRUE) {
  meta <- fread(metadata_path, sep = "\t", data.table = FALSE,
                na.strings = c("", "NA"))

  # Child-Pugh: string vazia é ausência, não B-C
  meta$child_pugh <- trimws(as.character(meta$child_pugh))
  meta$child_pugh_group <- ifelse(
    is.na(meta$child_pugh) | meta$child_pugh == "",
    NA_character_,
    ifelse(meta$child_pugh == "A", "A", "B-C")
  )
  meta$child_pugh_group <- factor(meta$child_pugh_group, levels = c("A", "B-C"))

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

  tss_counts <- table(meta$tss[meta$group != "GTEx_Normal"])
  tss_keep <- names(tss_counts)[tss_counts >= 10]
  meta$tss_collapsed <- factor(
    ifelse(meta$group == "GTEx_Normal", "GTEx",
           ifelse(meta$tss %in% tss_keep, as.character(meta$tss), "Outros")),
    levels = c(sort(tss_keep), "Outros", "GTEx")
  )

  meta$duplicata_gtex <- duplicated(meta$subject) & meta$group == "GTEx_Normal"
  n_duplicata <- sum(meta$duplicata_gtex)
  if (dedup_gtex) meta <- meta[!meta$duplicata_gtex, , drop = FALSE]

  for (v in names(meta)) if (is.factor(meta[[v]])) meta[[v]] <- droplevels(meta[[v]])

  list(
    meta = meta,
    n_duplicata = n_duplicata,
    tumor = meta[meta$group == "TCGA_Tumor", , drop = FALSE],
    adjacent = meta[meta$group == "TCGA_Adjacent", , drop = FALSE],
    gtex = meta[meta$group == "GTEx_Normal", , drop = FALSE],
    tcga = meta[meta$group %in% c("TCGA_Tumor", "TCGA_Adjacent"), , drop = FALSE],
    tcga_gtex = meta[meta$group %in% c("TCGA_Tumor", "GTEx_Normal"), , drop = FALSE]
  )
}

carregar_expressao <- function(expr_path) {
  expr_dt <- fread(expr_path, sep = "\t", data.table = FALSE)
  rownames(expr_dt) <- expr_dt$gene
  X <- as.matrix(expr_dt[, setdiff(names(expr_dt), "gene"), drop = FALSE])
  storage.mode(X) <- "numeric"
  X
}

log2tpm_para_tpm <- function(X) {
  tpm <- 2^X - 0.001
  tpm[tpm < 0] <- 0
  tpm
}

# Modelo limma com grupo de interesse + covariáveis. Retorna tabela completa.
rodar_limma <- function(X, d, group_col, levels, covariates = character(),
                        fdr_cut = 0.05, lfc_cut = 1, min_n = 6) {
  if (is.numeric(d[[group_col]])) {
    dd <- d[!is.na(d[[group_col]]), , drop = FALSE]
  } else {
    dd <- d[d[[group_col]] %in% levels, , drop = FALSE]
    dd[[group_col]] <- factor(as.character(dd[[group_col]]), levels = levels)
  }
  keep <- complete.cases(dd[, c(group_col, covariates), drop = FALSE])
  dd <- dd[keep, , drop = FALSE]
  for (v in names(dd)) if (is.factor(dd[[v]])) dd[[v]] <- droplevels(dd[[v]])
  if (nrow(dd) < min_n) return(NULL)
  fml <- as.formula(paste("~ 0 +", group_col,
                          if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""))
  design <- model.matrix(fml, data = dd)
  colnames(design) <- sub(group_col, "group", colnames(design), fixed = TRUE)
  if (qr(design)$rank < ncol(design)) return(NULL)
  fit <- lmFit(X[, dd$sample, drop = FALSE], design)
  cvec <- setNames(rep(0, ncol(design)), colnames(design))
  if (is.numeric(dd[[group_col]])) {
    coef_num <- c(paste0("group", group_col), "group")
    coef_num <- coef_num[coef_num %in% colnames(design)][1]
    if (is.na(coef_num)) return(NULL)
    cvec[coef_num] <- 1
  } else {
    if (!all(paste0("group", levels) %in% colnames(design))) return(NULL)
    cvec[paste0("group", levels[1])] <- -1
    cvec[paste0("group", levels[2])] <-  1
  }
  fit2 <- eBayes(contrasts.fit(fit, cvec), trend = TRUE, robust = TRUE)
  tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt$regulation <- ifelse(tt$adj.P.Val < fdr_cut & tt$logFC >  lfc_cut, "Up",
                   ifelse(tt$adj.P.Val < fdr_cut & tt$logFC < -lfc_cut, "Down", "NS"))
  tt
}

metricas_concordancia <- function(antes, depois) {
  if (is.null(antes) || is.null(depois)) {
    return(data.frame(degs_antes = NA, degs_depois = NA, degs_comuns = NA,
                      jaccard = NA, spearman_logfc = NA,
                      genes_antes_up = NA, genes_antes_down = NA,
                      genes_depois_up = NA, genes_depois_down = NA))
  }
  sa <- antes$gene[antes$regulation != "NS"]
  sb <- depois$gene[depois$regulation != "NS"]
  comuns <- intersect(antes$gene, depois$gene)
  data.frame(
    degs_antes = length(sa), degs_depois = length(sb),
    degs_comuns = length(intersect(sa, sb)),
    jaccard = if (length(union(sa, sb))) length(intersect(sa, sb)) / length(union(sa, sb)) else NA_real_,
    spearman_logfc = suppressWarnings(cor(
      antes$logFC[match(comuns, antes$gene)],
      depois$logFC[match(comuns, depois$gene)],
      method = "spearman", use = "complete.obs")),
    genes_antes_up = sum(antes$regulation == "Up"),
    genes_antes_down = sum(antes$regulation == "Down"),
    genes_depois_up = sum(depois$regulation == "Up"),
    genes_depois_down = sum(depois$regulation == "Down")
  )
}
