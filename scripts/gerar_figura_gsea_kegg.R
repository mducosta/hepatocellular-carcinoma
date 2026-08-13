suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
d <- read.csv("results/enrichment/GSEA_KEGG.csv", stringsAsFactors = FALSE)
d <- d[!is.na(d$ID) & d$ID != "" & !is.na(d$NES), ]
d <- d[!is.na(d$p.adjust) & d$p.adjust < 0.05, ]
cat("GSEA KEGG significativos:", nrow(d), "\n")
if (nrow(d) > 0) {
  d <- d %>% arrange(NES)
  d$Description <- gsub(" - ", " - ", d$Description)
  p <- ggplot(d, aes(x = reorder(Description, NES), y = NES, fill = NES < 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#4575B4", "FALSE" = "#D73027"),
                      guide = "none") +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(x = NULL, y = "Normalized Enrichment Score (NES)",
         title = "GSEA KEGG — vias significativamente enriquecidas (FDR < 0,05)",
         subtitle = "Genes da via hsa05417 ordenados por logFC | NES < 0 = down em LIHC") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
  ggsave("results/enrichment/GSEA_KEGG_barplot.png", p,
         width = 9, height = 4.5, dpi = 300)
  cat("Figura salva: results/enrichment/GSEA_KEGG_barplot.png\n")
}
