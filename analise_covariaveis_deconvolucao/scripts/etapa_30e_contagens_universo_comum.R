#!/usr/bin/env Rscript
# DESeq2 no mesmo universo gênico do TPM (12.401 genes), para comparação justa.
suppressPackageStartupMessages({library(data.table); library(DESeq2)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_out <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "sensibilidades")
dir_de  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "expressao_diferencial_covariaveis")
FDR <- 0.05; LFC <- 1

cnt <- fread(file.path(ROOT, "dados", "processed", "gdc_counts", "gdc_counts_proteincoding.tsv.gz"), data.table = FALSE)
mat <- as.matrix(cnt[, -1, drop = FALSE]); rownames(mat) <- cnt[[1]]
ann <- fread(file.path(ROOT, "dados", "processed", "gdc_counts", "gdc_counts_gene_annotation.tsv"), data.table = FALSE)
tpm <- fread(file.path(ROOT, "dados", "processed", "expressao_log2tpm_filtrada.tsv.gz"), data.table = FALSE)
genes_tpm <- tpm[[1]]
ann <- ann[ann$gene_name %in% genes_tpm, ]
ann <- ann[ann$gene_id %in% rownames(mat), ]
# escolhe, por símbolo, o gene_id de maior expressão média
media <- rowMeans(mat[ann$gene_id, , drop = FALSE])
ann$media <- media[match(ann$gene_id, names(media))]
ann <- ann[order(-ann$media), ]
ann <- ann[!duplicated(ann$gene_name), ]
mat <- mat[ann$gene_id, , drop = FALSE]
cat("universo comum:", nrow(mat), "genes\n")

meta <- fread(file.path(ROOT, "dados", "processed", "metadata_analise.tsv"), data.table = FALSE, na.strings = c("", "NA"))
meta <- meta[meta$sample %in% colnames(mat), , drop = FALSE]
meta <- meta[match(colnames(mat), meta$sample), , drop = FALSE]
meta$group <- factor(meta$group, levels = c("TCGA_Adjacent", "TCGA_Tumor"))
meta$sex <- factor(meta$sex); meta$age <- as.numeric(meta$age)
meta$stage_group <- factor(ifelse(meta$stage %in% c("I", "II"), "I-II",
                           ifelse(meta$stage %in% c("III", "IV"), "III-IV", NA)), levels = c("I-II", "III-IV"))
meta$grade_group <- factor(ifelse(meta$grade %in% c("G1", "G2"), "G1-G2",
                           ifelse(meta$grade %in% c("G3", "G4"), "G3-G4", NA)), levels = c("G1-G2", "G3-G4"))

rodar <- function(d, formula, contrast, nome) {
  d <- d[complete.cases(d[, all.vars(formula), drop = FALSE]), , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = mat[, d$sample, drop = FALSE], colData = d, design = formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = contrast, alpha = FDR)
  df <- as.data.frame(res); df$gene_id <- rownames(df)
  df <- merge(df, ann[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)
  df$gene <- ifelse(is.na(df$gene_name) | df$gene_name == "", df$gene_id, df$gene_name)
  df$regulation <- ifelse(df$padj < FDR & df$log2FoldChange > LFC, "Up",
                   ifelse(df$padj < FDR & df$log2FoldChange < -LFC, "Down", "NS"))
  fwrite(df, file.path(dir_out, paste0("contagens_comum_DESeq2_", nome, ".csv")))
  data.frame(metodo = "DESeq2_universo_comum", contraste = nome, n = nrow(d), genes = nrow(mat),
             degs = sum(df$regulation != "NS", na.rm = TRUE),
             up = sum(df$regulation == "Up", na.rm = TRUE),
             down = sum(df$regulation == "Down", na.rm = TRUE))
}
res <- rbind(
  rodar(meta, ~ group + sex + age, c("group", "TCGA_Tumor", "TCGA_Adjacent"), "tumor_vs_adjacente"),
  rodar(meta[meta$group == "TCGA_Tumor" & !is.na(meta$stage_group), ], ~ stage_group + sex + age,
        c("stage_group", "III-IV", "I-II"), "estagio"),
  rodar(meta[meta$group == "TCGA_Tumor" & !is.na(meta$grade_group), ], ~ grade_group + sex + age,
        c("grade_group", "G3-G4", "G1-G2"), "grau"))
fwrite(res, file.path(dir_out, "contagens_comum_resumo.csv"))
print(res)

pares <- list(tumor_vs_adjacente = "modelo_tumor_vs_adjacente_base_completo.csv",
              estagio = "modelo_estagio_base_completo.csv", grau = "modelo_grau_base_completo.csv")
conc <- list()
for (nm in names(pares)) {
  lim <- fread(file.path(dir_de, pares[[nm]]), data.table = FALSE)
  df <- fread(file.path(dir_out, paste0("contagens_comum_DESeq2_", nm, ".csv")), data.table = FALSE)
  df <- df[order(-df$baseMean), ]; df <- df[!duplicated(df$gene), ]
  com <- intersect(lim$gene, df$gene)
  a <- lim$gene[lim$regulation != "NS"]; b <- df$gene[df$regulation != "NS"]
  ov <- intersect(a, b)
  conc[[nm]] <- data.frame(contraste = nm, metodo = "DESeq2_universo_comum",
    genes_comuns = length(com), degs_limma = length(a), degs_metodo = length(b), degs_comuns = length(ov),
    jaccard = length(ov)/length(union(a, b)),
    spearman_logfc = cor(lim$logFC[match(com, lim$gene)], df$log2FoldChange[match(com, df$gene)],
                         method = "spearman", use = "complete.obs"),
    direcao_concordante = mean(sign(lim$logFC[match(ov, lim$gene)]) == sign(df$log2FoldChange[match(ov, df$gene)])))
}
conc <- do.call(rbind, conc)
fwrite(conc, file.path(dir_out, "contagens_comum_concordancia.csv"))
print(conc)
