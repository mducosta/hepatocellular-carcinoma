# =============================================================================
# Nanotec 2026 — GINAB/UENF
# Etapa 02 — Expressão diferencial (limma) + volcano + heatmap + correlação MMP9
# -----------------------------------------------------------------------------
# Os dados do painel JÁ estão em log2(TPM + 0.001) (auditado na Etapa 01:
# min = -9.97 = log2(0.001)). Portanto NÃO reaplicamos log2; usamos limma direto
# (lmFit + eBayes), adequado para dados contínuos em escala log2.
#
# CONTRASTES:
#   1. Tumor (TCGA-LIHC) vs Normal (GTEx Liver)  -> contraste principal
#      (ATENÇÃO: há confundimento TCGA/GTEx de plataforma aqui; por isso o 2.)
#   2. Tumor vs Adjacent (ambos TCGA)            -> validação sem lote de origem
#   3. Adjacent vs Normal                        -> contexto (TCGA vs GTEx)
#
# SAÍDAS (outputs/):
#   tables/  DEG_LIHC_vs_Normal.csv, DEG_LIHC_vs_Adjacent.csv, DEG_resumo.csv
#   volcano/ Volcano_LIHC_vs_Normal.png
#   heatmap/ Heatmap_painel_41genes.png
#   tables/  correlacao_MMP9_com_painel.csv
#
# Executar a partir da pasta nanotec2026/:
#   Rscript scripts/02_de_analysis.R
# =============================================================================

cat("============================================================\n")
cat("NANOTEC 2026 - Etapa 02: Expressão diferencial + figuras\n")
cat("============================================================\n")

base     <- "."
dir_proc <- file.path(base, "dados", "processed")
dir_tab  <- file.path(base, "outputs", "tables")
dir_vol  <- file.path(base, "outputs", "volcano")
dir_hm   <- file.path(base, "outputs", "heatmap")
for (d in c(dir_tab, dir_vol, dir_hm)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# pacotes
for (p in c("limma","ggplot2","ggrepel","pheatmap")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("limma")) { BiocManager::install(p, ask = FALSE, update = FALSE)
    } else install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pheatmap))

# -----------------------------------------------------------------------------
# 1. Leitura
# -----------------------------------------------------------------------------
mat <- read.delim(file.path(dir_proc, "liver_panel_41genes.tsv"),
                  header = TRUE, sep = "\t", check.names = FALSE, quote = "")
rownames(mat) <- mat[[1]]; mat <- mat[, -1, drop = FALSE]   # genes x amostras
annot <- read.delim(file.path(dir_proc, "samples_annotation.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE, quote = "")
m <- match(colnames(mat), annot[["sample"]])
stopifnot(!anyNA(m))
grupo <- annot[["grupo"]][m]
grupo <- factor(grupo, levels = c("Normal", "Adjacent", "Tumor"))
cat(sprintf("[DADOS] %d genes x %d amostras.\n", nrow(mat), ncol(mat)))
print(table(grupo))

# matriz para limma: genes x amostras (valores já em log2)
expr <- as.matrix(mat)
storage.mode(expr) <- "double"

# -----------------------------------------------------------------------------
# 2. DE com limma (3 grupos + contrastes)
# -----------------------------------------------------------------------------
cat("\n[limma] Ajustando modelo com 3 grupos e contrastes...\n")
design <- model.matrix(~ 0 + grupo)
colnames(design) <- levels(grupo)
fit <- lmFit(expr, design)
contr <- makeContrasts(
  Tumor_vs_Normal   = Tumor - Normal,
  Tumor_vs_Adjacent = Tumor - Adjacent,
  Adjacent_vs_Normal = Adjacent - Normal,
  levels = design
)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, trend = FALSE, robust = TRUE)

FC_CUT <- 1.0      # |log2FC| > 1
FDR_CUT <- 0.05    # FDR < 0.05

extract_de <- function(coef_name) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "none")
  tt <- data.frame(gene = rownames(tt), tt, stringsAsFactors = FALSE)
  tt[["signif"]] <- (!is.na(tt[["adj.P.Val"]])) & (tt[["adj.P.Val"]] < FDR_CUT) &
                    (abs(tt[["logFC"]]) > FC_CUT)
  tt[["direcao"]] <- ifelse(tt[["signif"]], ifelse(tt[["logFC"]] > 0, "Up", "Down"), "NS")
  tt <- tt[order(tt[["adj.P.Val"]]), ]
  tt
}

de_tumor_normal   <- extract_de("Tumor_vs_Normal")
de_tumor_adjacent <- extract_de("Tumor_vs_Adjacent")
de_adj_normal     <- extract_de("Adjacent_vs_Normal")

cat(sprintf("[DE] Tumor vs Normal  : %d DEGs (%d Up / %d Down)\n",
            sum(de_tumor_normal[["signif"]]),
            sum(de_tumor_normal[["signif"]] & de_tumor_normal[["logFC"]] > 0),
            sum(de_tumor_normal[["signif"]] & de_tumor_normal[["logFC"]] < 0)))
cat(sprintf("[DE] Tumor vs Adjacent: %d DEGs (%d Up / %d Down)\n",
            sum(de_tumor_adjacent[["signif"]]),
            sum(de_tumor_adjacent[["signif"]] & de_tumor_adjacent[["logFC"]] > 0),
            sum(de_tumor_adjacent[["signif"]] & de_tumor_adjacent[["logFC"]] < 0)))

write.csv(de_tumor_normal,   file.path(dir_tab, "DEG_LIHC_vs_Normal.csv"),   row.names = FALSE)
write.csv(de_tumor_adjacent, file.path(dir_tab, "DEG_LIHC_vs_Adjacent.csv"), row.names = FALSE)
write.csv(de_adj_normal,     file.path(dir_tab, "DEG_Adjacent_vs_Normal.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Volcano plot (Tumor vs Normal) com destaque MMP9 e MMPs
# -----------------------------------------------------------------------------
cat("\n[FIGURA] Volcano plot (Tumor vs Normal)...\n")
vd <- de_tumor_normal
vd[["nlog10fdr"]] <- -log10(vd[["adj.P.Val"]])
vd[["nlog10fdr"]][is.infinite(vd[["nlog10fdr"]])] <- max(vd[["nlog10fdr"]][is.finite(vd[["nlog10fdr"]])])

rotular <- vd[["gene"]] == "MMP9" |
           (vd[["signif"]] & grepl("^MMP", vd[["gene"]])) |
           (vd[["signif"]] & grepl("^TIMP", vd[["gene"]]))

p <- ggplot(vd, aes(x = logFC, y = nlog10fdr)) +
  geom_point(aes(color = direcao), alpha = 0.75, size = 2) +
  scale_color_manual(values = c("NS" = "grey70", "Up" = "#d73027", "Down" = "#4575b4"),
                     labels = c("NS", "Down", "Up"), name = "LIHC vs Normal") +
  geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = vd[rotular, ], aes(label = gene),
                  size = 3.4, max.overlaps = 25, min.segment.length = 0) +
  labs(title = "Volcano — painel ECM/MMP/TIMP (LIHC vs Normal)",
       x = "log2 fold-change", y = "-log10(FDR)") +
  theme_minimal(base_size = 12)
ggsave(file.path(dir_vol, "Volcano_LIHC_vs_Normal.png"), p, width = 8, height = 6, dpi = 150)

# -----------------------------------------------------------------------------
# 4. Heatmap dos 41 genes (z-score por gene)
# -----------------------------------------------------------------------------
cat("[FIGURA] Heatmap do painel...\n")
z <- t(scale(t(as.matrix(mat))))   # z-score por gene
z[!is.finite(z)] <- 0
# limita para visualização
z[z > 3] <- 3; z[z < -3] <- -3

# ordena amostras por grupo
ord <- order(as.integer(grupo))
annot_df <- data.frame(grupo = grupo, row.names = colnames(mat))

pheatmap(z[, ord], cluster_cols = FALSE, cluster_rows = TRUE,
         show_colnames = FALSE, fontsize_row = 8,
         annotation_col = annot_df,
         main = "Painel ECM/MMP/TIMP (41 genes) — z-score por gene",
         color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
         border_color = NA,
         filename = file.path(dir_hm, "Heatmap_painel_41genes.png"),
         width = 10, height = 7)

# -----------------------------------------------------------------------------
# 5. Correlação de MMP9 com os demais genes do painel
# -----------------------------------------------------------------------------
cat("\n[CORRELACAO] MMP9 vs demais genes (Spearman, amostras tumorais)...\n")
tumor_idx <- grupo == "Tumor"
m9 <- as.numeric(mat["MMP9", tumor_idx])
out_cor <- data.frame(gene = character(), rho = numeric(), p = numeric(),
                      stringsAsFactors = FALSE)
for (g in rownames(mat)) {
  if (g == "MMP9") next
  ct <- cor.test(m9, as.numeric(mat[g, tumor_idx]), method = "spearman",
                 exact = FALSE)
  out_cor <- rbind(out_cor, data.frame(gene = g, rho = ct$estimate, p = ct$p.value,
                                       stringsAsFactors = FALSE))
}
out_cor[["padj"]] <- p.adjust(out_cor[["p"]], method = "BH")
out_cor <- out_cor[order(out_cor[["padj"]]), ]
write.csv(out_cor, file.path(dir_tab, "correlacao_MMP9_com_painel.csv"), row.names = FALSE)
cat("[CORRELACAO] Top correlações com MMP9 (|rho|>0.3):\n")
print(head(out_cor[abs(out_cor[["rho"]]) > 0.3, c("gene","rho","padj")], 20))

# -----------------------------------------------------------------------------
# 6. Resumo consolidado
# -----------------------------------------------------------------------------
resumo <- data.frame(
  contraste = c("Tumor_vs_Normal","Tumor_vs_Adjacent","Adjacent_vs_Normal"),
  n_DEG = c(sum(de_tumor_normal[["signif"]]), sum(de_tumor_adjacent[["signif"]]),
            sum(de_adj_normal[["signif"]])),
  n_Up = c(sum(de_tumor_normal[["signif"]] & de_tumor_normal[["logFC"]]>0),
           sum(de_tumor_adjacent[["signif"]] & de_tumor_adjacent[["logFC"]]>0),
           sum(de_adj_normal[["signif"]] & de_adj_normal[["logFC"]]>0)),
  n_Down = c(sum(de_tumor_normal[["signif"]] & de_tumor_normal[["logFC"]]<0),
             sum(de_tumor_adjacent[["signif"]] & de_tumor_adjacent[["logFC"]]<0),
             sum(de_adj_normal[["signif"]] & de_adj_normal[["logFC"]]<0)),
  stringsAsFactors = FALSE
)
write.csv(resumo, file.path(dir_tab, "DEG_resumo.csv"), row.names = FALSE)

cat("\n[FIM] Etapa 02 concluída. Resultados em outputs/tables/, outputs/volcano/, outputs/heatmap/.\n")
