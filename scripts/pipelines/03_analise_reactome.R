# ===================================================================
# ENRIQUECIMENTO REACTOME — ORA + GSEA
# LIHC vs Normal — Via Lipid and Atherosclerosis (hsa05417)
# ===================================================================
# Entrada: outputs/deg/DEG_LIHC_vs_Normal_full.csv
# Saídas: outputs/enrichment/Reactome_*
# ===================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(rio)
})

PROJECT_ROOT <- normalizePath(".", mustWork = TRUE)
setwd(PROJECT_ROOT)

DEG_FILE <- "outputs/deg/DEG_LIHC_vs_Normal_full.csv"
stopifnot(file.exists(DEG_FILE))
deg <- read.csv(DEG_FILE, stringsAsFactors = FALSE)

map_df <- bitr(deg$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db, drop = TRUE)
deg <- deg %>% left_join(map_df, by = c("gene_symbol" = "SYMBOL"))

up_entrez   <- deg %>% filter(regulation == "Up_LIHC")   %>% pull(ENTREZID) %>% unique()
down_entrez <- deg %>% filter(regulation == "Down_LIHC") %>% pull(ENTREZID) %>% unique()

# ------------------------------------------------------------------
# 1) ORA — Reactome (genoma completo)
# ------------------------------------------------------------------
ora_reactome_up <- tryCatch(
  enrichPathway(gene = up_entrez, organism = "human",
                pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
  error = function(e) NULL)
ora_reactome_down <- tryCatch(
  enrichPathway(gene = down_entrez, organism = "human",
                pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
  error = function(e) NULL)

export_reactome <- function(obj, path) {
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    rio::export(as.data.frame(obj), path)
    cat(sprintf("  -> %s (%d termos)\n", path, nrow(as.data.frame(obj))))
  } else {
    rio::export(data.frame(), path)
    cat(sprintf("  -> %s (vazio)\n", path))
  }
}

cat("Reactome ORA up:\n");   export_reactome(ora_reactome_up,   "outputs/enrichment/Reactome_ORA_up.csv")
cat("Reactome ORA down:\n"); export_reactome(ora_reactome_down, "outputs/enrichment/Reactome_ORA_down.csv")

# ------------------------------------------------------------------
# 2) GSEA — Reactome (ranked by logFC)
# ------------------------------------------------------------------
ranked <- deg %>%
  filter(!is.na(ENTREZID), !is.na(logFC)) %>%
  arrange(desc(logFC)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gsea_genes <- ranked$logFC
names(gsea_genes) <- ranked$ENTREZID

gsea_reactome <- tryCatch(
  gsePathway(geneList = gsea_genes, organism = "human",
             eps = 0, pvalueCutoff = 0.1, seed = 4721),
  error = function(e) NULL)

if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
  gsea_df <- as.data.frame(gsea_reactome)
  gsea_df <- gsea_df[!is.na(gsea_df$ID), , drop = FALSE]
  rio::export(gsea_df, "outputs/enrichment/Reactome_GSEA.csv")
  cat(sprintf("Reactome GSEA: %d conjuntos (signif. padj<0.05: %d)\n",
              nrow(gsea_df), sum(gsea_df$p.adjust < 0.05, na.rm = TRUE)))
} else {
  rio::export(data.frame(), "outputs/enrichment/Reactome_GSEA.csv")
  cat("Reactome GSEA: vazio\n")
}

# ------------------------------------------------------------------
# 3) Figura — dotplot ORA up/down
# ------------------------------------------------------------------
plots <- list()
for (nm in c("up", "down")) {
  obj <- get(paste0("ora_reactome_", nm))
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    d <- as.data.frame(obj) %>% head(12) %>% mutate(set = paste0(nm, "_LIHC"))
    plots[[nm]] <- d
  }
}
if (length(plots) > 0) {
  rdf <- bind_rows(plots) %>% mutate(logFDR = -log10(p.adjust))
  p <- ggplot(rdf, aes(x = logFDR, y = reorder(Description, logFDR),
                       size = Count, color = logFDR)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ set, scales = "free_y", ncol = 1) +
    scale_color_gradient(low = "#74add1", high = "#d73027", name = "-log10(FDR)") +
    scale_size(range = c(3, 9), name = "Gene Count") +
    labs(x = "-log10(FDR)", y = "Reactome pathway",
         title = "Enriquecimento Reactome — ORA (up e down)") +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 8))
  ggsave("outputs/enrichment/Reactome_dotplot.png", p, width = 12, height = 9, dpi = 300)
  cat("Figura: outputs/enrichment/Reactome_dotplot.png\n")
}

cat("\n=== ENRIQUECIMENTO REACTOME CONCLUÍDO ===\n")
