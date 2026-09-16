# =============================================================================
# Nanotec 2026 — GINAB/UENF  (I Workshop Interdisciplinar de Nanotecnologia e Inovação)
# Etapa 05 — ANÁLISES ROBUSTAS COMPLEMENTARES (transcriptoma completo)
# -----------------------------------------------------------------------------
# Objetivo: ampliar a robustez da linha "remodelamento da MEC + MMP9 no HCC"
# com análises adicionais sobre o transcriptoma COMPLETO (36k genes):
#
#   1. DE genoma-wide (limma) Tumor vs Normal / Tumor vs Adjacent
#   2. GSEA (fgsea): Hallmark, KEGG, GO:BP  (assinatura global do TME)
#   3. ssGSEA de tipos celulares do fígado (atlas Aizarani, Nature 2019) e
#      de células imunes (marcadores canônicos) -> Tumor vs Normal + corr. c/ MMP9
#   4. Correlação de MMP9 com checkpoints imunes e marcadores de angiogênese
#   5. Sobrevida multi-desfecho (OS, DSS, PFI, DFI) + Cox multivariado (MMP9 + sexo)
#   6. Escore prognóstico ECM-MMP (MMP9+MMP1+MMP12+MMP14) -> KM + Cox
#   7. Estratificação por sexo (expressão de MMP9 + sobrevida)
#
# Saídas: outputs/robust/
#
# Executar a partir da pasta nanotec2026/:
#   Rscript scripts/05_analises_robustas.R
# =============================================================================

cat("============================================================\n")
cat("NANOTEC 2026 - Etapa 05: Análises robustas complementares\n")
cat("============================================================\n")

base     <- "."
dir_proc <- file.path(base, "dados", "processed")
dir_raw  <- file.path(base, "dados", "raw")
dir_out  <- file.path(base, "outputs", "robust")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(ggrepel)
  library(pheatmap); library(survival); library(survminer)
  library(fgsea); library(msigdbr); library(GSVA)
})
set.seed(2026)

# -----------------------------------------------------------------------------
# 0. Leitura dos dados
# -----------------------------------------------------------------------------
cat("\n[0] Lendo dados...\n")
t0 <- Sys.time()

# transcriptoma completo (símbolos), genes x amostras
full <- fread(file.path(dir_proc, "liver_transcriptome_symbols.tsv"),
              header = TRUE, sep = "\t", data.table = FALSE)
rownames(full) <- full[[1]]; full <- full[, -1, drop = FALSE]
cat(sprintf("  transcriptoma completo: %d genes x %d amostras\n", nrow(full), ncol(full)))

# painel 41 genes
panel <- read.delim(file.path(dir_proc, "liver_panel_41genes.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE, quote = "")
rownames(panel) <- panel[[1]]; panel <- panel[, -1, drop = FALSE]

# anotação
annot <- read.delim(file.path(dir_proc, "samples_annotation.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE, quote = "")
m <- match(colnames(full), annot[["sample"]])
stopifnot(!anyNA(m))
grupo <- factor(annot[["grupo"]][m], levels = c("Normal", "Adjacent", "Tumor"))
cat("  grupos:\n"); print(table(grupo))

# sexo (do fenótipo)
pheno <- fread(file.path(dir_raw, "TcgaTargetGTEX_phenotype.txt.gz"),
               header = TRUE, sep = "\t", data.table = FALSE, quote = "")
mp <- match(colnames(full), pheno[["sample"]])
sexo <- pheno[["_gender"]][mp]
sexo[!sexo %in% c("Female", "Male")] <- NA
cat(sprintf("  sexo: %d F / %d M / %d NA\n",
            sum(sexo == "Female", na.rm = TRUE),
            sum(sexo == "Male", na.rm = TRUE),
            sum(is.na(sexo))))

# sobrevivência
f_surv <- list.files(dir_raw, pattern = "TCGA_survival_data", full.names = TRUE)[1]
surv <- read.delim(f_surv, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")

expr_full <- as.matrix(full); storage.mode(expr_full) <- "double"

# -----------------------------------------------------------------------------
# 1. DE genoma-wide (limma)
# -----------------------------------------------------------------------------
cat("\n[1] DE genoma-wide (limma)...\n")
design <- model.matrix(~ 0 + grupo)
colnames(design) <- levels(grupo)
fit <- lmFit(expr_full, design)
contr <- makeContrasts(Tumor_vs_Normal = Tumor - Normal,
                       Tumor_vs_Adjacent = Tumor - Adjacent, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, trend = FALSE, robust = TRUE)

de_gw <- topTable(fit2, coef = "Tumor_vs_Normal", number = Inf, sort.by = "none")
de_gw <- data.frame(gene = rownames(de_gw), de_gw, stringsAsFactors = FALSE)
de_gw[["signif"]] <- (!is.na(de_gw[["adj.P.Val"]])) & de_gw[["adj.P.Val"]] < 0.05 &
                     abs(de_gw[["logFC"]]) > 1
cat(sprintf("  Tumor vs Normal: %d DEGs genoma-wide (%d up / %d down)\n",
            sum(de_gw[["signif"]]),
            sum(de_gw[["signif"]] & de_gw[["logFC"]] > 0),
            sum(de_gw[["signif"]] & de_gw[["logFC"]] < 0)))
de_gw <- de_gw[order(de_gw[["adj.P.Val"]]), ]
write.csv(de_gw, file.path(dir_out, "DEG_genomewide_Tumor_vs_Normal.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. GSEA (fgsea) — Hallmark, KEGG, GO:BP
# -----------------------------------------------------------------------------
cat("\n[2] GSEA (fgsea) — Hallmark, KEGG, GO:BP...\n")

rank_genes <- de_gw[["logFC"]]
names(rank_genes) <- de_gw[["gene"]]
rank_genes <- sort(rank_genes[is.finite(rank_genes)], decreasing = TRUE)

make_sets <- function(collection, subcollection = NULL) {
  d <- if (is.null(subcollection)) {
    msigdbr(species = "Homo sapiens", collection = collection)
  } else {
    msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
  }
  sets <- split(d$gene_symbol, d$gs_name)
  sets
}
run_fgsea <- function(sets, label) {
  res <- fgsea(pathways = sets, stats = rank_genes, minSize = 15, maxSize = 500,
               nPermSimple = 10000)
  res <- res[order(res$pval), ]
  res$leadingEdge <- vapply(res$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  as.data.frame(res)
}

# Hallmark
setsH <- make_sets("H")
fH <- run_fgsea(setsH, "Hallmark")
write.csv(fH, file.path(dir_out, "GSEA_Hallmark.csv"), row.names = FALSE)
cat(sprintf("  Hallmark: %d conjuntos, %d FDR<0.05\n", nrow(fH), sum(fH$padj < 0.05, na.rm=TRUE)))

# KEGG
setsK <- make_sets("C2", "CP:KEGG_LEGACY")
fK <- run_fgsea(setsK, "KEGG")
write.csv(fK, file.path(dir_out, "GSEA_KEGG.csv"), row.names = FALSE)
cat(sprintf("  KEGG: %d conjuntos, %d FDR<0.05\n", nrow(fK), sum(fK$padj < 0.05, na.rm=TRUE)))

# GO:BP
setsGO <- make_sets("C5", "GO:BP")
fGO <- run_fgsea(setsGO, "GO:BP")
write.csv(fGO, file.path(dir_out, "GSEA_GO_BP.csv"), row.names = FALSE)
cat(sprintf("  GO:BP: %d conjuntos, %d FDR<0.05\n", nrow(fGO), sum(fGO$padj < 0.05, na.rm=TRUE)))

# Figura: top vias Hallmark+KEGG (NES)
plot_gsea <- function(df, titulo, n = 24) {
  df <- df[!is.na(df$padj) & df$padj < 0.05 & !is.na(df$NES), , drop = FALSE]
  df <- df[!duplicated(df$pathway), , drop = FALSE]
  df$pathway <- as.character(df$pathway)
  df <- df[order(df$NES), , drop = FALSE]
  if (nrow(df) > n) {
    keep <- c(seq_len(n %/% 2), (nrow(df) - n %/% 2 + 1):nrow(df))
    df <- df[unique(keep), , drop = FALSE]
  }
  df$pathway <- factor(df$pathway, levels = df$pathway)
  ggplot(df, aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"), guide = "none") +
    geom_vline(xintercept = 0, color = "black") +
    labs(title = titulo, x = "Normalized Enrichment Score (NES)", y = "") +
    theme_minimal(base_size = 11)
}
ggsave(file.path(dir_out, "GSEA_Hallmark_top.png"),
       plot_gsea(fH, "Hallmark — top vias (LIHC vs Normal)"),
       width = 9, height = 7, dpi = 150)
ggsave(file.path(dir_out, "GSEA_KEGG_top.png"),
       plot_gsea(fK, "KEGG — top vias (LIHC vs Normal)"),
       width = 9, height = 7, dpi = 150)

# -----------------------------------------------------------------------------
# 3. Deconvolução de tipos celulares (ssGSEA)
# -----------------------------------------------------------------------------
cat("\n[3] ssGSEA — tipos celulares do fígado (Aizarani) + imunes...\n")

# (a) Atlas hepático Aizarani (Nature 2019) — agrupa clusters por tipo celular
c8 <- msigdbr(species = "Homo sapiens", collection = "C8")
az <- c8[grepl("^AIZARANI_LIVER_", c8$gs_name), ]
az$celltype <- sub("^AIZARANI_LIVER_C[0-9]+_", "", az$gs_name)
az$celltype <- sub("_[0-9]+$", "", az$celltype)
az$celltype <- gsub("_", " ", az$celltype)
sets_az <- split(az$gene_symbol, az$celltype)
cat(sprintf("  Aizarani: %d tipos celulares hepáticos\n", length(sets_az)))

# (b) Células imunes (marcadores canônicos, Bindea-like)
immune_markers <- list(
  "CD8 T"        = c("CD8A", "CD8B"),
  "CD4 T"        = c("CD4"),
  "Treg"         = c("FOXP3", "IL2RA", "CTLA4"),
  "NK"           = c("NKG7", "KLRD1", "KLRK1", "NCR1"),
  "B cells"      = c("CD19", "CD79A", "MS4A1"),
  "Plasma cells" = c("SDC1", "JCHAIN", "MZB1"),
  "Macrofago M1" = c("NOS2", "IL1B", "TNF", "CXCL10", "IRF5"),
  "Macrofago M2" = c("CD163", "MRC1", "MSR1", "ARG1"),
  "Monocito"     = c("CD14", "FCGR3A", "LYZ"),
  "Neutrofilo"   = c("CEACAM8", "FCGR3B", "S100A8", "S100A9", "CSF3R"),
  "Dendritica"   = c("ITGAX", "CD1C", "CLEC9A", "LILRA4", "CD83"),
  "Mastocito"    = c("KIT", "CPA3", "TPSAB1", "MS4A2"),
  "Endotelio"    = c("PECAM1", "CDH5", "VWF", "ENG"),
  "Fibroblasto/CAF" = c("ACTA2", "PDGFRB", "FAP", "DCN", "LUM")
)

sets_all <- c(sets_az, immune_markers)
# mantém só conjuntos com genes presentes
sets_all <- lapply(sets_all, function(g) intersect(g, rownames(expr_full)))
sets_all <- sets_all[lengths(sets_all) >= 3]
cat(sprintf("  Total de conjuntos para ssGSEA: %d\n", length(sets_all)))

param <- ssgseaParam(expr_full, sets_all)
ss <- gsva(param)
ss <- as.data.frame(ss)
write.csv(data.frame(celltype = rownames(ss), ss, check.names = FALSE),
          file.path(dir_out, "ssGSEA_celltype_scores.csv"), row.names = FALSE)

# Tumor vs Normal (Wilcoxon) por tipo
tumor_idx <- grupo == "Tumor"; normal_idx <- grupo == "Normal"
cmp <- data.frame(celltype = rownames(ss), p = NA_real_, diff = NA_real_,
                  median_Tumor = NA_real_, median_Normal = NA_real_, stringsAsFactors = FALSE)
for (i in seq_len(nrow(ss))) {
  a <- as.numeric(ss[i, tumor_idx]); b <- as.numeric(ss[i, normal_idx])
  w <- suppressWarnings(wilcox.test(a, b))
  cmp$p[i] <- w$p.value
  cmp$diff[i] <- median(a) - median(b)
  cmp$median_Tumor[i] <- median(a); cmp$median_Normal[i] <- median(b)
}
cmp$padj <- p.adjust(cmp$p, method = "BH")
cmp <- cmp[order(cmp$padj), ]
write.csv(cmp, file.path(dir_out, "ssGSEA_celltype_Tumor_vs_Normal.csv"), row.names = FALSE)
cat("  Top tipos com diferença Tumor vs Normal (FDR<0.05):\n")
print(head(cmp[cmp$padj < 0.05, c("celltype","diff","median_Tumor","median_Normal","padj")], 15))

# Correlação MMP9 (tumor) com cada score celular
m9_t <- as.numeric(panel["MMP9", colnames(panel)[grupo == "Tumor"]])
cor_cell <- data.frame(celltype = rownames(ss), rho = NA_real_, p = NA_real_, stringsAsFactors = FALSE)
for (i in seq_len(nrow(ss))) {
  ct <- cor.test(m9_t, as.numeric(ss[i, tumor_idx]), method = "spearman", exact = FALSE)
  cor_cell$rho[i] <- ct$estimate; cor_cell$p[i] <- ct$p.value
}
cor_cell$padj <- p.adjust(cor_cell$p, method = "BH")
cor_cell <- cor_cell[order(-cor_cell$rho), ]
write.csv(cor_cell, file.path(dir_out, "correlacao_MMP9_celltype.csv"), row.names = FALSE)
cat("\n  Correlação de MMP9 com tipos celulares (tumor, Spearman):\n")
print(head(cor_cell[cor_cell$padj < 0.05, c("celltype","rho","padj")], 15))

# Heatmap dos scores celulares (top 30 por significância, amostras ordenadas)
top_cells <- head(cmp$celltype[cmp$padj < 0.05], 30)
if (length(top_cells) >= 2) {
  z <- t(scale(t(as.matrix(ss[top_cells, , drop = FALSE]))))
  z[!is.finite(z)] <- 0; z[z > 3] <- 3; z[z < -3] <- -3
  ord <- order(as.integer(grupo))
  annot_df <- data.frame(grupo = grupo, row.names = colnames(expr_full))
  pheatmap(z[, ord], cluster_cols = FALSE, cluster_rows = TRUE,
           show_colnames = FALSE, fontsize_row = 8, annotation_col = annot_df,
           main = "Scores celulares (ssGSEA) — Tumor vs Normal",
           color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
           border_color = NA,
           filename = file.path(dir_out, "Heatmap_celltype_scores.png"),
           width = 10, height = 9)
}

# -----------------------------------------------------------------------------
# 4. Correlação MMP9 com checkpoints imunes e angiogênese
# -----------------------------------------------------------------------------
cat("\n[4] MMP9 x checkpoints imunes / angiogênese...\n")
checkpoints <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
                 "IDO1","ENTPD1","NT5E","VSIR","CD276","LGALS9","SIGLEC15",
                 "CD80","CD86","CD40","ICOS","TNFRSF9","TNFRSF4","TNFRSF18")
angiog <- c("VEGFA","KDR","FLT1","PECAM1","ANGPT2","TEK","NOTCH1","DLL4","EFNB2")
marcadores <- unique(c(checkpoints, angiog))
marcadores <- marcadores[marcadores %in% rownames(expr_full)]

cor_ck <- data.frame(gene = marcadores, rho = NA_real_, p = NA_real_, stringsAsFactors = FALSE)
for (g in marcadores) {
  ct <- cor.test(m9_t, as.numeric(expr_full[g, colnames(full)[grupo == "Tumor"]]),
                 method = "spearman", exact = FALSE)
  cor_ck$rho[match(g, cor_ck$gene)] <- ct$estimate
  cor_ck$p[match(g, cor_ck$gene)] <- ct$p.value
}
cor_ck$padj <- p.adjust(cor_ck$p, method = "BH")
cor_ck$classe <- ifelse(cor_ck$gene %in% checkpoints, "Checkpoint imune", "Angiogênese")
cor_ck <- cor_ck[order(-cor_ck$rho), ]
write.csv(cor_ck, file.path(dir_out, "correlacao_MMP9_checkpoints_angiogenese.csv"), row.names = FALSE)
cat("  Correlações com FDR<0.05:\n")
print(cor_ck[cor_ck$padj < 0.05, c("gene","classe","rho","padj")])

# -----------------------------------------------------------------------------
# 5. Sobrevida multi-desfecho + Cox multivariado (MMP9 + sexo)
# -----------------------------------------------------------------------------
cat("\n[5] Sobrevida multi-desfecho + Cox multivariado...\n")
tumor_samp <- colnames(full)[grupo == "Tumor"]
idx <- match(tumor_samp, surv[["sample"]])

endpoints <- c("OS", "DSS", "PFI", "DFI")
genes_alvo <- c("MMP9", "MMP1", "MMP12", "MMP14")

res_surv <- data.frame()
for (endp in endpoints) {
  tcol <- paste0(endp, ".time"); ecol <- endp
  ok <- !is.na(idx) & !is.na(surv[[tcol]][idx]) & !is.na(surv[[ecol]][idx])
  time <- surv[[tcol]][idx[ok]]; evt <- surv[[ecol]][idx[ok]]
  # sexo das amostras tumorais
  sx <- sexo[match(tumor_samp, colnames(full))][ok]
  for (g in genes_alvo) {
    x <- as.numeric(expr_full[g, tumor_samp][ok])
    # univariado
    c1 <- coxph(Surv(time, evt) ~ x)
    # multivariado (x + sexo)
    dfmv <- data.frame(x = x, sx = factor(sx, levels = c("Male","Female")))
    c2 <- tryCatch(coxph(Surv(time, evt) ~ x + sx, data = dfmv),
                   error = function(e) NULL)
    res_surv <- rbind(res_surv, data.frame(
      endpoint = endp, gene = g,
      HR_univ = exp(coef(c1))[1],
      p_univ = summary(c1)$coefficients[1, "Pr(>|z|)"],
      HR_multi = if (is.null(c2)) NA else exp(coef(c2))["x"],
      p_multi = if (is.null(c2)) NA else summary(c2)$coefficients["x", "Pr(>|z|)"],
      stringsAsFactors = FALSE))
  }
}
write.csv(res_surv, file.path(dir_out, "sobrevida_multidesfecho_Cox.csv"), row.names = FALSE)
cat("  Resultados (univariado, HR por unidade de log2):\n")
print(res_surv[res_surv$gene == "MMP9", c("endpoint","HR_univ","p_univ","HR_multi","p_multi")])

# -----------------------------------------------------------------------------
# 6. Escore prognóstico ECM-MMP (z-score médio MMP9+MMP1+MMP12+MMP14)
# -----------------------------------------------------------------------------
cat("\n[6] Escore prognóstico ECM-MMP...\n")
gs <- intersect(genes_alvo, rownames(expr_full))
zmat <- t(scale(t(expr_full[gs, tumor_samp])))
escore <- colMeans(zmat, na.rm = TRUE)

ok <- !is.na(idx) & !is.na(surv[["OS.time"]][idx]) & !is.na(surv[["OS"]][idx])
time <- surv[["OS.time"]][idx[ok]]; evt <- surv[["OS"]][idx[ok]]
es <- escore[ok]
grp <- factor(ifelse(es >= median(es, na.rm = TRUE), "High", "Low"), levels = c("Low","High"))
sfit <- survfit(Surv(time, evt) ~ grp)
lr <- survdiff(Surv(time, evt) ~ grp); p_lr <- 1 - pchisq(lr$chisq, df = 1)
cox_es <- coxph(Surv(time, evt) ~ es)
png(file.path(dir_out, "KM_ECM_MMP_score_OS.png"), width = 900, height = 700, res = 130)
plt <- ggsurvplot(sfit, data = data.frame(time, evt, grp), pval = TRUE, risk.table = TRUE,
                  title = "Sobrevida global — escore ECM-MMP (MMP9+MMP1+MMP12+MMP14)",
                  xlab = "Tempo (dias)", ylab = "Probabilidade de sobrevida",
                  legend.title = "Escore", palette = c("#4575b4", "#d73027"))
print(plt); dev.off()
cat(sprintf("  Escore ECM-MMP (OS): log-rank p=%.4g, Cox HR=%.3f p=%.4g\n",
            p_lr, exp(coef(cox_es))[1], summary(cox_es)$coefficients[1, "Pr(>|z|)"]))

# -----------------------------------------------------------------------------
# 7. Estratificação por sexo
# -----------------------------------------------------------------------------
cat("\n[7] Estratificação por sexo...\n")
sx_t <- sexo[match(tumor_samp, colnames(full))]
m9_all <- as.numeric(expr_full["MMP9", tumor_samp])
sx_t <- factor(sx_t, levels = c("Female", "Male"))
res_sex <- data.frame()
for (s in c("Female", "Male")) {
  keep <- which(sx_t == s)
  other <- which(sx_t != s)
  w <- wilcox.test(m9_all[keep], m9_all[other])
  # sobrevida MMP9 dentro do sexo
  ok2 <- !is.na(idx) & !is.na(surv[["OS.time"]][idx]) & !is.na(surv[["OS"]][idx])
  k <- keep[keep %in% which(ok2)]
  x2 <- as.numeric(expr_full["MMP9", tumor_samp][k])
  t2 <- surv[["OS.time"]][idx[k]]; e2 <- surv[["OS"]][idx[k]]
  c2 <- tryCatch(coxph(Surv(t2, e2) ~ x2), error = function(e) NULL)
  res_sex <- rbind(res_sex, data.frame(sexo = s,
    n = length(keep),
    mediana_MMP9 = median(m9_all[keep]),
    p_expressao = w$p.value,
    HR_MMP9 = if (is.null(c2)) NA else exp(coef(c2))[1],
    p_MMP9 = if (is.null(c2)) NA else summary(c2)$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE))
}
write.csv(res_sex, file.path(dir_out, "estratificacao_sexo_MMP9.csv"), row.names = FALSE)
print(res_sex)

# boxplot MMP9 por grupo x sexo
dfp <- data.frame(MMP9 = as.numeric(expr_full["MMP9", ]),
                  grupo = grupo,
                  sexo = factor(sexo, levels = c("Female","Male")))
dfp <- dfp[!is.na(dfp$sexo), ]
p <- ggplot(dfp, aes(x = grupo, y = MMP9, fill = sexo)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  scale_fill_manual(values = c("Female" = "#d95f02", "Male" = "#7570b3")) +
  labs(title = "MMP9 por grupo e sexo", x = "", y = "log2(TPM+0.001)") +
  theme_minimal(base_size = 12)
ggsave(file.path(dir_out, "Boxplot_MMP9_grupo_sexo.png"), p, width = 7, height = 5, dpi = 150)

cat(sprintf("\n[FIM] Etapa 05 concluída em %.1f min. Resultados em outputs/robust/.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
