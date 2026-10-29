#!/usr/bin/env Rscript
# Concordância das contagens (DESeq2/edgeR) com os modelos limma.
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_out <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "sensibilidades")
dir_de  <- file.path(ROOT, "resultados", "analises_avancadas", "tabelas", "expressao_diferencial_covariaveis")
pares <- list(tumor_vs_adjacente = "modelo_tumor_vs_adjacente_base_completo.csv",
              estagio = "modelo_estagio_base_completo.csv",
              grau = "modelo_grau_base_completo.csv")
conc <- list()
for (nm in names(pares)) {
  lim <- fread(file.path(dir_de, pares[[nm]]), data.table = FALSE)
  for (met in c("DESeq2", "edgeR")) {
    f <- file.path(dir_out, paste0("contagens_", met, "_", nm, ".csv"))
    if (!file.exists(f)) next
    df <- fread(f, data.table = FALSE)
    ord <- if ("baseMean" %in% names(df)) df$baseMean else df$logCPM
    df <- df[order(-ord), ]; df <- df[!duplicated(df$gene), ]
    com <- intersect(lim$gene, df$gene)
    a <- lim$gene[lim$regulation != "NS"]; b <- df$gene[df$regulation != "NS"]
    ov <- intersect(a, b)
    conc[[paste(nm, met)]] <- data.frame(
      contraste = nm, metodo = met, genes_comuns = length(com),
      degs_limma = length(a), degs_metodo = length(b), degs_comuns = length(ov),
      jaccard = length(ov) / length(union(a, b)),
      spearman_logfc = cor(lim$logFC[match(com, lim$gene)], df$log2FoldChange[match(com, df$gene)],
                           method = "spearman", use = "complete.obs"),
      direcao_concordante = mean(sign(lim$logFC[match(ov, lim$gene)]) ==
                                 sign(df$log2FoldChange[match(ov, df$gene)])))
  }
}
conc <- do.call(rbind, conc)
fwrite(conc, file.path(dir_out, "contagens_concordancia_limma.csv"))
print(conc)
