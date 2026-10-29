#!/usr/bin/env Rscript
# =============================================================================
# etapa_30b_contagens_DESeq2_edgeR.R
# Sensibilidade com contagens brutas STAR (GDC): DESeq2 e edgeR.
# Compara os DEGs com os modelos limma em TPM.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(DESeq2); library(edgeR)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_cnt <- file.path(ROOT, "dados", "processed", "gdc_counts")
dir_out <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "sensibilidades")
dir_de  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "expressao_diferencial_covariaveis")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
FDR <- 0.05; LFC <- 1

cnt <- fread(file.path(dir_cnt, "gdc_counts_proteincoding.tsv.gz"), data.table = FALSE)
mat <- as.matrix(cnt[, -1, drop = FALSE]); rownames(mat) <- cnt[[1]]
ann <- fread(file.path(dir_cnt, "gdc_counts_gene_annotation.tsv"), data.table = FALSE)
meta <- fread(file.path(ROOT, "dados", "processed", "metadata_analise.tsv"), data.table = FALSE,
              na.strings = c("", "NA"))
meta <- meta[meta$sample %in% colnames(mat), , drop = FALSE]
meta <- meta[match(colnames(mat), meta$sample), , drop = FALSE]
meta$group <- factor(meta$group, levels = c("TCGA_Adjacent", "TCGA_Tumor"))
meta$sex <- factor(meta$sex, levels = c("female", "male"))
meta$age <- as.numeric(meta$age)
meta$stage_group <- factor(ifelse(meta$stage %in% c("I", "II"), "I-II",
                           ifelse(meta$stage %in% c("III", "IV"), "III-IV", NA)),
                           levels = c("I-II", "III-IV"))
meta$grade_group <- factor(ifelse(meta$grade %in% c("G1", "G2"), "G1-G2",
                           ifelse(meta$grade %in% c("G3", "G4"), "G3-G4", NA)),
                           levels = c("G1-G2", "G3-G4"))
cat("amostras:", ncol(mat), "genes:", nrow(mat), "\n")

add_ann <- function(df) {
  df$gene_id <- rownames(df)
  df <- merge(df, ann, by = "gene_id", all.x = TRUE)
  df$gene <- ifelse(is.na(df$gene_name) | df$gene_name == "", df$gene_id, df$gene_name)
  df
}
classificar <- function(df, lfc_col = "log2FoldChange", p_col = "padj") {
  df$regulation <- ifelse(df[[p_col]] < FDR & df[[lfc_col]] > LFC, "Up",
                   ifelse(df[[p_col]] < FDR & df[[lfc_col]] < -LFC, "Down", "NS"))
  df
}

rodar_deseq <- function(d, formula, contrast, nome) {
  d <- d[complete.cases(d[, all.vars(formula), drop = FALSE]), , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = mat[, d$sample, drop = FALSE],
                                colData = d, design = formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = contrast, alpha = FDR, independentFiltering = TRUE)
  df <- as.data.frame(res); df <- add_ann(df); df <- classificar(df)
  fwrite(df, file.path(dir_out, paste0("contagens_DESeq2_", nome, ".csv")))
  data.frame(metodo = "DESeq2", contraste = nome, n = nrow(d),
             degs = sum(df$regulation != "NS"), up = sum(df$regulation == "Up"),
             down = sum(df$regulation == "Down"))
}

rodar_edger <- function(d, formula, coef, nome) {
  d <- d[complete.cases(d[, all.vars(formula), drop = FALSE]), , drop = FALSE]
  y <- DGEList(counts = mat[, d$sample, drop = FALSE], group = d[[all.vars(formula)[1]]])
  y <- calcNormFactors(y)
  design <- model.matrix(formula, data = d)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = coef)
  tt <- topTags(qlf, n = Inf)$table
  tt$log2FoldChange <- tt$logFC / log(2)
  tt$padj <- tt$FDR
  df <- add_ann(tt); df <- classificar(df)
  fwrite(df, file.path(dir_out, paste0("contagens_edgeR_", nome, ".csv")))
  data.frame(metodo = "edgeR", contraste = nome, n = nrow(d),
             degs = sum(df$regulation != "NS"), up = sum(df$regulation == "Up"),
             down = sum(df$regulation == "Down"))
}

res <- list()
# Tumor vs adjacente
res[[1]] <- rodar_deseq(meta, ~ group + sex + age, c("group", "TCGA_Tumor", "TCGA_Adjacent"), "tumor_vs_adjacente")
res[[2]] <- rodar_edger(meta, ~ group + sex + age, "groupTCGA_Tumor", "tumor_vs_adjacente")

# Estágio (tumores)
d_st <- meta[meta$group == "TCGA_Tumor" & !is.na(meta$stage_group), , drop = FALSE]
res[[3]] <- rodar_deseq(d_st, ~ stage_group + sex + age, c("stage_group", "III-IV", "I-II"), "estagio")
res[[4]] <- rodar_edger(d_st, ~ stage_group + sex + age, "stage_groupIII-IV", "estagio")

# Grau (tumores)
d_gr <- meta[meta$group == "TCGA_Tumor" & !is.na(meta$grade_group), , drop = FALSE]
res[[5]] <- rodar_deseq(d_gr, ~ grade_group + sex + age, c("grade_group", "G3-G4", "G1-G2"), "grau")
res[[6]] <- rodar_edger(d_gr, ~ grade_group + sex + age, "grade_groupG3-G4", "grau")

resumo <- do.call(rbind, res)
fwrite(resumo, file.path(dir_out, "contagens_resumo.csv"))
print(resumo)

# Concordância com limma
pares <- list(
  tumor_vs_adjacente = "modelo_tumor_vs_adjacente_base_completo.csv",
  estagio = "modelo_estagio_base_completo.csv",
  grau = "modelo_grau_base_completo.csv"
)
conc <- list()
for (nm in names(pares)) {
  lim <- fread(file.path(dir_de, pares[[nm]]), data.table = FALSE)
  for (met in c("DESeq2", "edgeR")) {
    df <- fread(file.path(dir_out, paste0("contagens_", met, "_", nm, ".csv")), data.table = FALSE)
    # colapsa duplicatas de símbolo pela maior expressão média
    ord <- if ("baseMean" %in% names(df)) df$baseMean else df$logCPM
    df <- df[order(-ord), ]; df <- df[!duplicated(df$gene), ]
    com <- intersect(lim$gene, df$gene)
    a <- lim$gene[lim$regulation != "NS"]; b <- df$gene[df$regulation != "NS"]
    conc[[paste(nm, met)]] <- data.frame(
      contraste = nm, metodo = met, genes_comuns = length(com),
      degs_limma = length(a), degs_metodo = length(b), degs_comuns = length(intersect(a, b)),
      jaccard = length(intersect(a, b)) / length(union(a, b)),
      spearman_logfc = cor(lim$logFC[match(com, lim$gene)], df$log2FoldChange[match(com, df$gene)],
                           method = "spearman", use = "complete.obs"),
      direcao_concordante = mean(sign(lim$logFC[match(intersect(a, b), lim$gene)]) ==
                                 sign(df$log2FoldChange[match(intersect(a, b), df$gene)])))
  }
}
conc <- do.call(rbind, conc)
fwrite(conc, file.path(dir_out, "contagens_concordancia_limma.csv"))
print(conc)
