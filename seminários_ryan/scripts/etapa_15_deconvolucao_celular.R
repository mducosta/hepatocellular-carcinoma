# =============================================================================
# etapa_15_deconvolucao_celular.R
# -----------------------------------------------------------------------------
# Deconvolução celular de transcriptomas bulk de fígado.
#
# Três métodos complementares:
#   1. EPIC (frações; inclui otherCells, que em fígado representa sobretudo
#      hepatócitos, além de CAFs, endotélio, macrófagos e linfócitos).
#   2. quanTIseq (frações de 10 tipos imunes + Other).
#   3. xCell (escores de enriquecimento de 64 tipos celulares, incluindo
#      Hepatocytes, Endothelial cells, Fibroblasts e subtipos imunes).
#
# Entradas: log2(TPM + 0,001). Para EPIC/quanTIseq, converte-se para TPM
# linear (2^x - 0,001).
#
# Saídas: frações/escores, QC, comparação entre grupos, associações clínicas,
# concordância entre métodos, Cox ajustado e figuras.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(EPIC)
  library(quantiseqr)
  library(xCell)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))

dir_proc <- file.path(ROOT, "dados", "processed")
dir_out  <- file.path(ROOT, "resultados", "analises_avancadas")
dir_tab  <- file.path(dir_out, "tabelas", "deconvolucao")
dir_fig  <- file.path(dir_out, "figuras", "deconvolucao")
dir_qc   <- file.path(dir_out, "qc")
for (d in c(dir_tab, dir_fig, dir_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

prep <- preparar_metadados(file.path(dir_proc, "metadata_analise.tsv"), dedup_gtex = TRUE)
meta <- prep$meta
X <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))
X <- X[, meta$sample, drop = FALSE]
tpm <- log2tpm_para_tpm(X)

log_lines <- c("DECONVOLUÇÃO CELULAR — LOG", "=============================",
               paste("Data:", Sys.time()), paste("R:", R.version.string),
               paste("Amostras:", ncol(X), "| Genes:", nrow(X)),
               paste("GTEx duplicado excluído:", prep$n_duplicata),
               "")

# -----------------------------------------------------------------------------
# 1. EPIC
# -----------------------------------------------------------------------------
message("Rodando EPIC...")
epic <- EPIC::EPIC(bulk = tpm, mRNA_cell = NULL, scaleExprs = TRUE,
                   withOtherCells = TRUE)
epic_frac <- as.data.frame(epic$cellFractions)
epic_frac$sample <- rownames(epic_frac)
rownames(epic_frac) <- NULL
epic_frac <- epic_frac[, c("sample", setdiff(names(epic_frac), "sample"))]
names(epic_frac)[-1] <- paste0("EPIC_", names(epic_frac)[-1])
fwrite(epic_frac, file.path(dir_tab, "fracoes_EPIC.csv"))

epic_gof <- epic$fit.gof
epic_gof$sample <- rownames(epic_gof)
rownames(epic_gof) <- NULL
fwrite(epic_gof, file.path(dir_tab, "QC_EPIC_convergencia.csv"))
log_lines <- c(log_lines,
               paste("EPIC: amostras =", nrow(epic_frac),
                     "| convergência código 0 =", sum(epic_gof$convergeCode == 0),
                     "| códigos != 0 =", sum(epic_gof$convergeCode != 0)),
               paste("EPIC: soma das frações (min/med/max) =",
                     paste(round(range(rowSums(epic_frac[, -1])), 3), collapse = " / ")))

# -----------------------------------------------------------------------------
# 2. quanTIseq
# -----------------------------------------------------------------------------
message("Rodando quanTIseq...")
qt <- quantiseqr::run_quantiseq(expression_data = tpm,
                               signature_matrix = "TIL10",
                               is_tumordata = TRUE)
qt <- as.data.frame(qt)
names(qt)[1] <- "sample"
names(qt)[-1] <- paste0("quanTIseq_", names(qt)[-1])
fwrite(qt, file.path(dir_tab, "fracoes_quanTIseq.csv"))
log_lines <- c(log_lines,
               paste("quanTIseq: amostras =", nrow(qt),
                     "| soma das frações (min/med/max) =",
                     paste(round(range(rowSums(qt[, -1])), 3), collapse = " / ")))

# -----------------------------------------------------------------------------
# 3. xCell
# -----------------------------------------------------------------------------
message("Rodando xCell...")
xc <- xCell::xCellAnalysis(X, rnaseq = TRUE, parallel.sz = 1)
xc_mat <- t(xc)
xc_df <- data.frame(sample = rownames(xc_mat), xc_mat, check.names = FALSE)
names(xc_df)[-1] <- paste0("xCell_", names(xc_df)[-1])
fwrite(xc_df, file.path(dir_tab, "escores_xCell.csv"))
log_lines <- c(log_lines,
               paste("xCell: amostras =", nrow(xc_df),
                     "| tipos/escores =", ncol(xc_df) - 1))

# -----------------------------------------------------------------------------
# 4. Tabela unificada
# -----------------------------------------------------------------------------
uni <- Reduce(function(a, b) merge(a, b, by = "sample", all = FALSE),
              list(meta[, c("sample", "subject", "group", "condition", "study",
                            "sex", "age", "stage", "stage_group", "grade",
                            "grade_group", "viral_status", "child_pugh_group",
                            "histology_group", "afp", "afp_group2",
                            "fibrosis_ishak", "fibrosis_adv2",
                            "vascular_invasion2", "tss", "tss_collapsed",
                            "OS", "OS_time", "vital_status")],
                   epic_frac, qt, xc_df))
fwrite(uni, file.path(dir_tab, "deconvolucao_unificada_amostras.csv"))

# -----------------------------------------------------------------------------
# 5. QC: correlação com marcadores conhecidos
# -----------------------------------------------------------------------------
marcadores <- list(
  Hepatocytes = c("ALB", "APOA1", "HP", "TTR", "CYP3A4", "CYP2E1"),
  Cholangiocytes = c("KRT19", "EPCAM", "KRT7", "SOX9"),
  Endothelial = c("PECAM1", "CDH5", "VWF", "KDR"),
  Fibroblasts = c("COL1A1", "COL1A2", "ACTA2", "PDGFRB"),
  Macrophages = c("CD68", "CD14", "CSF1R", "AIF1"),
  `T cells` = c("CD3D", "CD3E", "CD2", "TRAC"),
  `B cells` = c("CD19", "MS4A1", "CD79A"),
  NK = c("NKG7", "GNLY", "KLRD1"),
  Neutrophils = c("FCGR3B", "S100A8", "S100A9")
)

colunas_celulas <- c(
  "EPIC_otherCells", "xCell_Hepatocytes",
  "xCell_Cholangiocytes", "xCell_Epithelial cells",
  "EPIC_Endothelial", "xCell_Endothelial cells",
  "EPIC_CAFs", "xCell_Fibroblasts",
  "EPIC_Macrophages", "xCell_Macrophages", "quanTIseq_Macrophages.M1", "quanTIseq_Macrophages.M2",
  "EPIC_CD8_Tcells", "xCell_CD8+ T-cells", "quanTIseq_T.cells.CD8",
  "EPIC_CD4_Tcells", "xCell_CD4+ T-cells", "quanTIseq_T.cells.CD4",
  "EPIC_Bcells", "xCell_B-cells", "quanTIseq_B.cells",
  "EPIC_NKcells", "xCell_NK cells", "quanTIseq_NK.cells",
  "quanTIseq_Neutrophils", "xCell_Neutrophils"
)
colunas_celulas <- colunas_celulas[colunas_celulas %in% names(uni)]

cor_qc <- list()
for (ct in colunas_celulas) {
  for (mk in unique(unlist(marcadores))) {
    if (!mk %in% rownames(X)) next
    rho <- suppressWarnings(cor(uni[[ct]], X[mk, uni$sample],
                                method = "spearman", use = "complete.obs"))
    cor_qc[[paste(ct, mk)]] <- data.frame(
      celula = ct, marcador = mk, spearman = rho,
      familia = names(marcadores)[vapply(marcadores, function(z) mk %in% z, logical(1))][1],
      stringsAsFactors = FALSE)
  }
}
cor_qc <- do.call(rbind, cor_qc)
fwrite(cor_qc, file.path(dir_tab, "QC_correlacao_celulas_marcadores.csv"))

# Matriz resumida: melhor correlação esperada por célula
melhor_cor <- do.call(rbind, lapply(split(cor_qc, cor_qc$celula), function(d) {
  d[which.max(abs(d$spearman)), ]
}))
fwrite(melhor_cor, file.path(dir_tab, "QC_melhor_correlacao_por_celula.csv"))

# -----------------------------------------------------------------------------
# 6. Concordância entre métodos
# -----------------------------------------------------------------------------
pares_metodos <- data.frame(
  conceito = c("Hepatócitos/other", "Endotélio", "Fibroblastos/CAFs",
               "Macrófagos", "CD8 T", "CD4 T", "B", "NK"),
  metodo_a = c("EPIC_otherCells", "EPIC_Endothelial", "EPIC_CAFs",
               "EPIC_Macrophages", "EPIC_CD8_Tcells", "EPIC_CD4_Tcells",
               "EPIC_Bcells", "EPIC_NKcells"),
  metodo_b = c("xCell_Hepatocytes", "xCell_Endothelial cells", "xCell_Fibroblasts",
               "xCell_Macrophages", "xCell_CD8+ T-cells", "xCell_CD4+ T-cells",
               "xCell_B-cells", "xCell_NK cells"),
  metodo_c = c(NA, NA, NA, "quanTIseq_Macrophages.M1",
               "quanTIseq_T.cells.CD8", "quanTIseq_T.cells.CD4",
               "quanTIseq_B.cells", "quanTIseq_NK.cells"),
  stringsAsFactors = FALSE
)
conc <- list()
for (i in seq_len(nrow(pares_metodos))) {
  p <- pares_metodos[i, ]
  a <- p$metodo_a; b <- p$metodo_b; cc <- p$metodo_c
  if (!a %in% names(uni) || !b %in% names(uni)) next
  rho_ab <- suppressWarnings(cor(uni[[a]], uni[[b]], method = "spearman", use = "complete.obs"))
  rho_ac <- if (!is.na(cc) && cc %in% names(uni)) {
    suppressWarnings(cor(uni[[a]], uni[[cc]], method = "spearman", use = "complete.obs"))
  } else NA_real_
  rho_bc <- if (!is.na(cc) && cc %in% names(uni)) {
    suppressWarnings(cor(uni[[b]], uni[[cc]], method = "spearman", use = "complete.obs"))
  } else NA_real_
  conc[[i]] <- data.frame(conceito = p$conceito, metodo_a = a, metodo_b = b,
                          metodo_c = cc, spearman_ab = rho_ab,
                          spearman_ac = rho_ac, spearman_bc = rho_bc,
                          stringsAsFactors = FALSE)
}
conc <- do.call(rbind, conc)
fwrite(conc, file.path(dir_tab, "concordancia_metodos_deconvolucao.csv"))

# -----------------------------------------------------------------------------
# 7. Comparação entre grupos
# -----------------------------------------------------------------------------
grupos <- c("TCGA_Tumor", "TCGA_Adjacent", "GTEx_Normal")
comp_group <- list()
for (ct in colunas_celulas) {
  d <- uni[!is.na(uni[[ct]]) & !is.na(uni$group), ]
  if (nrow(d) < 20) next
  kw <- kruskal.test(d[[ct]] ~ d$group)
  med <- aggregate(d[[ct]], list(group = d$group), median, na.rm = TRUE)
  names(med)[2] <- "mediana"
  meanv <- aggregate(d[[ct]], list(group = d$group), mean, na.rm = TRUE)
  names(meanv)[2] <- "media"
  # pares não pareados
  p_tg <- tryCatch(wilcox.test(d[[ct]][d$group == "TCGA_Tumor"],
                               d[[ct]][d$group == "GTEx_Normal"])$p.value,
                   error = function(e) NA_real_)
  p_ag <- tryCatch(wilcox.test(d[[ct]][d$group == "TCGA_Adjacent"],
                               d[[ct]][d$group == "GTEx_Normal"])$p.value,
                   error = function(e) NA_real_)
  p_ta <- tryCatch(wilcox.test(d[[ct]][d$group == "TCGA_Tumor"],
                               d[[ct]][d$group == "TCGA_Adjacent"])$p.value,
                   error = function(e) NA_real_)
  # pareado tumor-adjacente
  pair_samples <- intersect(
    uni$subject[uni$group == "TCGA_Tumor"],
    uni$subject[uni$group == "TCGA_Adjacent"]
  )
  p_par <- NA_real_
  if (length(pair_samples) >= 10) {
    a <- uni[uni$group == "TCGA_Tumor" & uni$subject %in% pair_samples, c("subject", ct)]
    b <- uni[uni$group == "TCGA_Adjacent" & uni$subject %in% pair_samples, c("subject", ct)]
    m <- merge(a, b, by = "subject")
    if (nrow(m) >= 10) p_par <- tryCatch(wilcox.test(m[[2]], m[[3]], paired = TRUE)$p.value,
                                         error = function(e) NA_real_)
  }
  comp_group[[ct]] <- data.frame(
    celula = ct,
    mediana_tumor = med$mediana[med$group == "TCGA_Tumor"][1],
    mediana_adjacente = med$mediana[med$group == "TCGA_Adjacent"][1],
    mediana_gtex = med$mediana[med$group == "GTEx_Normal"][1],
    media_tumor = meanv$media[meanv$group == "TCGA_Tumor"][1],
    media_adjacente = meanv$media[meanv$group == "TCGA_Adjacent"][1],
    media_gtex = meanv$media[meanv$group == "GTEx_Normal"][1],
    p_kruskal = kw$p.value,
    p_tumor_vs_gtex = p_tg,
    p_adjacente_vs_gtex = p_ag,
    p_tumor_vs_adjacente = p_ta,
    p_pareado_tumor_adjacente = p_par,
    n = nrow(d), stringsAsFactors = FALSE
  )
}
comp_group <- do.call(rbind, comp_group)
comp_group$fdr_kruskal <- p.adjust(comp_group$p_kruskal, "BH")
comp_group$fdr_tumor_vs_gtex <- p.adjust(comp_group$p_tumor_vs_gtex, "BH")
comp_group$fdr_adjacente_vs_gtex <- p.adjust(comp_group$p_adjacente_vs_gtex, "BH")
comp_group$fdr_tumor_vs_adjacente <- p.adjust(comp_group$p_tumor_vs_adjacente, "BH")
comp_group$fdr_pareado <- p.adjust(comp_group$p_pareado_tumor_adjacente, "BH")
fwrite(comp_group, file.path(dir_tab, "comparacao_celulas_por_grupo.csv"))

# -----------------------------------------------------------------------------
# 8. Associações clínicas (somente tumor)
# -----------------------------------------------------------------------------
tum <- uni[uni$group == "TCGA_Tumor", , drop = FALSE]

assoc_var <- function(ct, var, tipo = c("cat", "num")) {
  tipo <- match.arg(tipo)
  d <- tum[!is.na(tum[[ct]]) & !is.na(tum[[var]]), , drop = FALSE]
  if (nrow(d) < 15) return(NULL)
  if (tipo == "num") {
    rho <- suppressWarnings(cor(d[[ct]], d[[var]], method = "spearman"))
    p <- suppressWarnings(cor.test(d[[ct]], d[[var]], method = "spearman")$p.value)
    med <- NA_real_
    n_grupos <- NA_integer_
    return(data.frame(celula = ct, variavel = var, tipo = tipo,
                      mediana = med, n_grupos = n_grupos,
                      spearman = rho, p = p, n = nrow(d), stringsAsFactors = FALSE))
  }
  g <- factor(as.character(d[[var]]))
  if (nlevels(g) < 2) return(NULL)
  p <- tryCatch(kruskal.test(d[[ct]] ~ g)$p.value, error = function(e) NA_real_)
  med <- tapply(d[[ct]], g, median, na.rm = TRUE)
  data.frame(celula = ct, variavel = var, tipo = tipo,
             mediana = paste(names(med), round(med, 4), sep = "=", collapse = "; "),
             n_grupos = nlevels(g), spearman = NA_real_, p = p, n = nrow(d),
             stringsAsFactors = FALSE)
}

vars_cat <- c("sex", "stage_group", "grade_group", "viral_status",
              "child_pugh_group", "histology_group", "afp_group2",
              "fibrosis_adv2", "vascular_invasion2")
vars_num <- c("age", "stage_num", "grade_num")

assoc <- list()
for (ct in colunas_celulas) {
  for (v in vars_cat) {
    if (!v %in% names(tum)) next
    a <- assoc_var(ct, v, "cat")
    if (!is.null(a)) assoc[[paste(ct, v)]] <- a
  }
  for (v in vars_num) {
    if (!v %in% names(tum)) next
    a <- assoc_var(ct, v, "num")
    if (!is.null(a)) assoc[[paste(ct, v)]] <- a
  }
}
assoc <- do.call(rbind, assoc)
if (!is.null(assoc)) {
  assoc$fdr_variavel <- ave(assoc$p, assoc$variavel, FUN = function(p) p.adjust(p, "BH"))
  assoc$fdr_global <- p.adjust(assoc$p, "BH")
  fwrite(assoc, file.path(dir_tab, "associacoes_celulas_clinica.csv"))
}

# -----------------------------------------------------------------------------
# 9. Sobrevida: Cox ajustado por idade, sexo, estágio e grau
# -----------------------------------------------------------------------------
cox_res <- list()
tum_os <- tum[!is.na(tum$OS) & !is.na(tum$OS_time) & tum$OS_time > 0, ]
for (ct in colunas_celulas) {
  d <- tum_os[!is.na(tum_os[[ct]]), ]
  d <- d[complete.cases(d[, c("OS", "OS_time", "age", "sex",
                              "stage_group", "grade_group")]), ]
  if (nrow(d) < 40 || sum(d$OS == 1) < 15) next
  d$score_z <- as.numeric(scale(d[[ct]]))
  fit <- tryCatch(coxph(Surv(OS_time, OS) ~ score_z + age + sex + stage_group + grade_group,
                        data = d), error = function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  co <- s$coefficients["score_z", ]
  ci <- s$conf.int["score_z", ]
  cox_res[[ct]] <- data.frame(
    celula = ct, HR_por_SD = ci[1], IC95_inf = ci[3], IC95_sup = ci[4],
    p = co["Pr(>|z|)"], n = fit$n, eventos = fit$nevent,
    stringsAsFactors = FALSE
  )
}
cox_res <- do.call(rbind, cox_res)
if (!is.null(cox_res)) {
  cox_res$fdr <- p.adjust(cox_res$p, "BH")
  fwrite(cox_res, file.path(dir_tab, "cox_celulas_ajustado.csv"))
}

# -----------------------------------------------------------------------------
# 10. Figuras
# -----------------------------------------------------------------------------
# 10a. EPIC por grupo
epic_long <- melt(as.data.table(uni[, c("group", grep("^EPIC_", names(uni), value = TRUE))]),
                  id.vars = "group", variable.name = "celula", value.name = "fracao")
epic_long <- epic_long[!is.na(epic_long$group), ]
p1 <- ggplot(epic_long, aes(x = group, y = fracao, fill = group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.65) +
  facet_wrap(~celula, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(TCGA_Tumor = "#D1495B", TCGA_Adjacent = "#EDAE49",
                               GTEx_Normal = "#00798C")) +
  labs(x = NULL, y = "Fração estimada (EPIC)", fill = "Grupo",
       title = "Deconvolução EPIC por grupo biológico") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(size = 8))
ggsave(file.path(dir_fig, "figura_fracoes_EPIC_por_grupo.png"), p1,
       width = 12, height = 7, dpi = 300)

# 10b. xCell key por grupo
xcell_key <- c("xCell_Hepatocytes", "xCell_Endothelial cells", "xCell_Fibroblasts",
               "xCell_Macrophages", "xCell_CD8+ T-cells", "xCell_NK cells",
               "xCell_B-cells", "xCell_Epithelial cells")
xcell_key <- xcell_key[xcell_key %in% names(uni)]
xc_long <- melt(as.data.table(uni[, c("group", xcell_key)]), id.vars = "group",
                variable.name = "celula", value.name = "escore")
xc_long <- xc_long[!is.na(xc_long$group), ]
p2 <- ggplot(xc_long, aes(x = group, y = escore, fill = group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.65) +
  facet_wrap(~celula, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(TCGA_Tumor = "#D1495B", TCGA_Adjacent = "#EDAE49",
                               GTEx_Normal = "#00798C")) +
  labs(x = NULL, y = "Escore xCell", fill = "Grupo",
       title = "Escores xCell por grupo biológico") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(size = 8))
ggsave(file.path(dir_fig, "figura_escores_xCell_por_grupo.png"), p2,
       width = 12, height = 7, dpi = 300)

# 10c. Hepatócitos por sexo e grupo
hep <- uni[!is.na(uni$sex) & !is.na(uni$group), ]
hep_long <- rbind(
  data.frame(group = hep$group, sex = hep$sex, valor = hep$xCell_Hepatocytes,
             metodo = "xCell Hepatocytes"),
  data.frame(group = hep$group, sex = hep$sex, valor = hep$EPIC_otherCells,
             metodo = "EPIC otherCells")
)
p3 <- ggplot(hep_long, aes(x = sex, y = valor, fill = sex)) +
  geom_boxplot(outlier.size = 0.4, width = 0.65) +
  facet_grid(metodo ~ group, scales = "free_y") +
  scale_fill_manual(values = c(female = "#7B2CBF", male = "#1B998B")) +
  labs(x = NULL, y = "Estimativa celular", fill = "Sexo",
       title = "Hepatócitos/otherCells por sexo e compartimento") +
  theme_bw(base_size = 10)
ggsave(file.path(dir_fig, "figura_hepatocitos_sexo.png"), p3,
       width = 9, height = 6, dpi = 300)

# 10d. Heatmap de concordância entre métodos
conc_mat <- as.matrix(conc[, c("spearman_ab", "spearman_ac", "spearman_bc")])
rownames(conc_mat) <- conc$conceito
pheatmap(conc_mat, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Concordância entre métodos (Spearman)",
         filename = file.path(dir_fig, "figura_concordancia_metodos.png"),
         width = 7, height = 5)

# 10e. Heatmap de associações clínicas (rho ou -log10 p)
if (!is.null(assoc)) {
  assoc_cat <- assoc[assoc$tipo == "cat", ]
  assoc_cat$logp <- -log10(assoc_cat$p)
  mat_assoc <- reshape(assoc_cat[, c("celula", "variavel", "logp")],
                       idvar = "celula", timevar = "variavel", direction = "wide")
  rownames(mat_assoc) <- mat_assoc$celula
  mat_assoc$celula <- NULL
  names(mat_assoc) <- sub("logp\\.", "", names(mat_assoc))
  mat_assoc[is.na(mat_assoc)] <- 0
  pheatmap(as.matrix(mat_assoc), cluster_rows = TRUE, cluster_cols = TRUE,
           color = colorRampPalette(c("white", "#FDBF6F", "#E31A1C"))(100),
           main = "Associação célula × clínica (-log10 P)",
           filename = file.path(dir_fig, "figura_associacoes_celulas_clinica.png"),
           width = 11, height = 8)
}

# 10f. Forest plot Cox
if (!is.null(cox_res)) {
  cox_res$celula <- factor(cox_res$celula, levels = cox_res$celula[order(cox_res$HR_por_SD)])
  p4 <- ggplot(cox_res, aes(x = HR_por_SD, y = celula)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_pointrange(aes(xmin = IC95_inf, xmax = IC95_sup,
                        colour = p < 0.05), size = 0.4) +
    scale_x_log10() +
    scale_colour_manual(values = c(`TRUE` = "#D1495B", `FALSE` = "grey50")) +
    labs(x = "HR por 1 SD (IC 95%, escala log)", y = NULL,
         title = "Cox ajustado: composição celular e sobrevida global",
         colour = "P < 0,05") +
    theme_bw(base_size = 10)
  ggsave(file.path(dir_fig, "figura_cox_celulas.png"), p4,
         width = 8, height = 7, dpi = 300)
}

# 10g. Kaplan-Meier para hepatócitos e macrófagos
km_key <- c("xCell_Hepatocytes", "EPIC_Macrophages", "xCell_CD8+ T-cells", "xCell_Fibroblasts")
for (ct in km_key) {
  if (!ct %in% names(tum_os)) next
  d <- tum_os[!is.na(tum_os[[ct]]), ]
  if (nrow(d) < 40) next
  d$grupo_expr <- factor(ifelse(d[[ct]] >= median(d[[ct]], na.rm = TRUE), "Alto", "Baixo"),
                         levels = c("Baixo", "Alto"))
  fit <- survfit(Surv(OS_time, OS) ~ grupo_expr, data = d)
  p <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = FALSE,
                  palette = c("#00798C", "#D1495B"),
                  legend.title = ct, legend.labs = c("Baixo", "Alto"),
                  xlab = "Tempo (dias)", ylab = "Sobrevida global",
                  title = paste("Sobrevida por", ct))
  ggsave(file.path(dir_fig, paste0("KM_", gsub("[^A-Za-z0-9]+", "_", ct), ".png")),
         p$plot, width = 6, height = 5, dpi = 300)
}

# -----------------------------------------------------------------------------
# 11. Resumo
# -----------------------------------------------------------------------------
top_celulas <- comp_group[order(comp_group$p_kruskal), c("celula", "p_kruskal", "fdr_kruskal")]
log_lines <- c(log_lines, "",
               "Top células por diferença entre grupos (Kruskal):",
               paste0("  ", top_celulas$celula[1:min(12, nrow(top_celulas))],
                      ": p=", signif(top_celulas$p_kruskal[1:min(12, nrow(top_celulas))], 3),
                      " | FDR=", signif(top_celulas$fdr_kruskal[1:min(12, nrow(top_celulas))], 3)),
               "",
               "Concordância entre métodos (Spearman):",
               paste0("  ", conc$conceito, ": EPIC-xCell=", round(conc$spearman_ab, 3),
                      ifelse(is.na(conc$spearman_ac), "",
                             paste0("; EPIC-quanTIseq=", round(conc$spearman_ac, 3),
                                    "; xCell-quanTIseq=", round(conc$spearman_bc, 3)))),
               "",
               "Cox ajustado (top 10 por p):",
               paste0("  ", cox_res$celula[order(cox_res$p)][1:min(10, nrow(cox_res))],
                      ": HR=", round(cox_res$HR_por_SD[order(cox_res$p)][1:min(10, nrow(cox_res))], 2),
                      " | p=", signif(cox_res$p[order(cox_res$p)][1:min(10, nrow(cox_res))], 3),
                      " | FDR=", signif(cox_res$fdr[order(cox_res$p)][1:min(10, nrow(cox_res))], 3)))
writeLines(log_lines, file.path(dir_qc, "log_deconvolucao_celular.txt"))

cat("Concluído.\n")
print(comp_group[order(comp_group$p_kruskal)[1:15],
                 c("celula", "mediana_tumor", "mediana_adjacente", "mediana_gtex",
                   "p_kruskal", "fdr_kruskal")])
