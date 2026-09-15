# =============================================================================
# Nanotec 2026 — GINAB/UENF
# Etapa 04 — Scores de enriquecimento por amostra (GSVA/ssGSEA) de vias ECM
# -----------------------------------------------------------------------------
# Usa o transcriptoma COMPLETO das amostras de fígado (símbolos HUGO, 36301 genes)
# para calcular escores de enriquecimento por amostra (ssGSEA via GSVA) com:
#   - Hallmark (50 conjuntos, incluindo EMT, Angiogenesis, Coagulation);
#   - KEGG focados em ECM (ECM-receptor interaction, Focal adhesion,
#     Protein digestion and absorption).
# Depois compara Tumor (LIHC) vs Normal (GTEx) por conjunto (teste t + BH).
#
# JUSTIFICATIVA: GSVA/ssGSEA exige o ranking de TODOS os genes da amostra como
# fundo — por isso usamos o transcriptoma completo (Etapa 00), não só o painel.
#
# Saídas (outputs/enrichment/):
#   GSVA_scores.csv, GSVA_summary.csv, GSVA_barplot.png
#
# Executar a partir da pasta nanotec2026/:
#   Rscript scripts/04_enrichment.R
# =============================================================================

cat("============================================================\n")
cat("NANOTEC 2026 - Etapa 04: GSVA/ssGSEA de vias ECM\n")
cat("============================================================\n")

base     <- "."
dir_proc <- file.path(base, "dados", "processed")
dir_enr  <- file.path(base, "outputs", "enrichment")
dir.create(dir_enr, showWarnings = FALSE, recursive = TRUE)

for (p in c("GSVA", "msigdbr", "data.table", "ggplot2")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "GSVA") BiocManager::install(p, ask = FALSE, update = FALSE)
    else install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# -----------------------------------------------------------------------------
# 1. Leitura do transcriptoma completo (símbolos)
# -----------------------------------------------------------------------------
cat("[DADOS] Lendo transcriptoma completo (símbolos)...\n")
dt <- fread(file.path(dir_proc, "liver_transcriptome_symbols.tsv"), header = TRUE)
setnames(dt, 1, "SYMBOL")
gene_sym <- dt[["SYMBOL"]]
samples  <- setdiff(names(dt), "SYMBOL")

# matriz genes x amostras
expr <- as.matrix(dt[, ..samples])
rownames(expr) <- gene_sym
storage.mode(expr) <- "double"
cat(sprintf("[DADOS] %d genes x %d amostras.\n", nrow(expr), ncol(expr)))

annot <- read.delim(file.path(dir_proc, "samples_annotation.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE, quote = "")
m <- match(samples, annot[["sample"]])
stopifnot(!anyNA(m))
grupo <- annot[["grupo"]][m]
grupo <- factor(grupo, levels = c("Normal", "Adjacent", "Tumor"))

# -----------------------------------------------------------------------------
# 2. Genesets: Hallmark + KEGG ECM
# -----------------------------------------------------------------------------
cat("[GENESETS] Carregando Hallmark e KEGG (msigdbr)...\n")
hall <- msigdbr(species = "Homo sapiens", collection = "H")
hall_list <- split(hall[["gene_symbol"]], hall[["gs_name"]])

kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
ecm_kegg <- c("KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_FOCAL_ADHESION")
kegg_list <- split(kegg[["gene_symbol"]], kegg[["gs_name"]])
kegg_list <- kegg_list[names(kegg_list) %in% ecm_kegg]

gs_all <- c(hall_list, kegg_list)
cat(sprintf("[GENESETS] %d conjuntos (Hallmark %d + KEGG ECM %d).\n",
            length(gs_all), length(hall_list), length(kegg_list)))

# -----------------------------------------------------------------------------
# 3. GSVA (ssGSEA)
# -----------------------------------------------------------------------------
cat("[GSVA] Calculando escores ssGSEA por amostra...\n")
param <- ssgseaParam(expr, gs_all, minSize = 10, maxSize = 500)
scores <- gsva(param)
# scores: genesets x amostras
cat(sprintf("[GSVA] %d conjuntos x %d amostras calculados.\n", nrow(scores), ncol(scores)))

write.csv(as.data.frame(scores), file.path(dir_enr, "GSVA_scores.csv"))

# -----------------------------------------------------------------------------
# 4. Tumor vs Normal por conjunto (teste t + BH)
# -----------------------------------------------------------------------------
idx_tumor <- which(grupo == "Tumor")
idx_norm  <- which(grupo == "Normal")

res <- data.frame(set = rownames(scores), stringsAsFactors = FALSE)
res[["media_Tumor"]] <- rowMeans(scores[, idx_tumor])
res[["media_Normal"]] <- rowMeans(scores[, idx_norm])
res[["diff"]] <- res[["media_Tumor"]] - res[["media_Normal"]]
res[["p"]] <- apply(scores, 1, function(r) {
  tryCatch(t.test(r[idx_tumor], r[idx_norm])$p.value, error = function(e) NA)
})
res[["padj"]] <- p.adjust(res[["p"]], method = "BH")
res <- res[order(res[["padj"]]), ]

write.csv(res, file.path(dir_enr, "GSVA_summary.csv"), row.names = FALSE)

cat("\n[GSVA] Conjuntos mais diferenciais (Tumor vs Normal, FDR<0.05):\n")
sig <- res[res[["padj"]] < 0.05 & !is.na(res[["padj"]]), ]
print(sig[, c("set", "diff", "padj")], row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Barplot dos top conjuntos
# -----------------------------------------------------------------------------
top <- head(res, 20)
top[["set"]][order(top[["diff"]])]
top[["set"]] <- factor(top[["set"]], levels = top[["set"]][order(top[["diff"]])])
p <- ggplot(top, aes(x = diff, y = set, fill = diff > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                    labels = c("Down em LIHC", "Up em LIHC"), name = "LIHC vs Normal") +
  labs(title = "GSVA/ssGSEA — top conjuntos (Tumor vs Normal)",
       x = "Diferença média de escore (Tumor - Normal)", y = "") +
  theme_minimal(base_size = 11)
ggsave(file.path(dir_enr, "GSVA_barplot.png"), p, width = 9, height = 6, dpi = 150)

cat("\n[FIM] Etapa 04 concluída. Resultados em outputs/enrichment/.\n")
