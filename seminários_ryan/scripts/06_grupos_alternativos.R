# =============================================================================
# 06_grupos_alternativos.R
# -----------------------------------------------------------------------------
# Inventário e exploração de grupos além de "tumor vs adjacente":
#   - tipos de amostra (incluindo recidivas);
#   - subtipos histológicos (HCC, mixed, fibrolamellar);
#   - grupos virais (HBV, HCV, HBV+HCV, none/other, unknown);
#   - centros TCGA (TSS);
#   - marcadores por subtipo histológico (descritivo; n pequeno).
#
# Saídas:
#   resultados/tabelas/21_inventario_grupos.csv
#   resultados/tabelas/22_histologia_marcadores.csv
#   resultados/figuras/26_grupos_dados.png
#   resultados/figuras/27_histologia_marcadores.png
#   resultados/qc/resumo_grupos_alternativos.txt
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_proc <- file.path(ROOT, "dados", "processed")
dir_tab  <- file.path(ROOT, "resultados", "tabelas")
dir_fig  <- file.path(ROOT, "resultados", "figuras")
dir_qc   <- file.path(ROOT, "resultados", "qc")
for (d in c(dir_tab, dir_fig, dir_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
              sep = "\t", data.table = FALSE, na.strings = c("", "NA"))
meta_all <- fread(file.path(dir_proc, "metadata_amostras.tsv"),
                  sep = "\t", data.table = FALSE, na.strings = c("", "NA"))
tcga <- meta[meta$group == "TCGA_Tumor", ]

# -----------------------------------------------------------------------------
# 1. Inventário de grupos
# -----------------------------------------------------------------------------
inv <- list()
add_inv <- function(categoria, grupo, n, total = NULL) {
  if (is.null(total)) total <- sum(n)
  inv[[length(inv) + 1]] <<- data.frame(
    categoria = categoria, grupo = grupo, n = n,
    percentual = round(100 * n / total, 1),
    stringsAsFactors = FALSE
  )
}

# Tipos de amostra (inclui recidivas)
tab_st <- as.data.frame(table(meta_all$sample_type), stringsAsFactors = FALSE)
names(tab_st) <- c("grupo", "n")
add_inv("Tipo de amostra", tab_st$grupo, tab_st$n, sum(tab_st$n))

# Grupos da análise principal
tab_gr <- as.data.frame(table(meta$group), stringsAsFactors = FALSE)
names(tab_gr) <- c("grupo", "n")
add_inv("Grupo analítico", tab_gr$grupo, tab_gr$n, sum(tab_gr$n))

# Histologia (TCGA tumor)
tab_hist <- as.data.frame(table(tcga$histology), stringsAsFactors = FALSE)
names(tab_hist) <- c("grupo", "n")
tab_hist <- tab_hist[tab_hist$n > 0, ]
add_inv("Histologia (TCGA tumor)", tab_hist$grupo, tab_hist$n, sum(tab_hist$n))

# Viral (TCGA tumor)
tab_vir <- as.data.frame(table(tcga$viral_status), stringsAsFactors = FALSE)
names(tab_vir) <- c("grupo", "n")
tab_vir <- tab_vir[tab_vir$n > 0, ]
add_inv("Sorologia viral (TCGA tumor)", tab_vir$grupo, tab_vir$n, sum(tab_vir$n))

# TSS (TCGA tumor)
tab_tss <- as.data.frame(table(tcga$tss), stringsAsFactors = FALSE)
names(tab_tss) <- c("grupo", "n")
tab_tss <- tab_tss[order(-tab_tss$n), ]
add_inv("Centro TCGA (TSS)", tab_tss$grupo, tab_tss$n, sum(tab_tss$n))

inventario <- do.call(rbind, inv)
write.csv(inventario, file.path(dir_tab, "21_inventario_grupos.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Figura de grupos
# -----------------------------------------------------------------------------
cores_tipo <- c(
  "Primary Tumor" = "#C0392B",
  "Solid Tissue Normal" = "#E67E22",
  "Normal Tissue" = "#2471A3",
  "Recurrent Tumor" = "#7D3C98"
)
pA <- ggplot(tab_st, aes(x = reorder(grupo, -n), y = n, fill = grupo)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.8) +
  scale_fill_manual(values = cores_tipo) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "A. Tipos de amostra",
       subtitle = "Inclui as 2 recidivas excluídas da análise principal",
       x = NULL, y = "n", fill = "Tipo") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid.minor = element_blank())

pB <- ggplot(tab_hist, aes(x = reorder(grupo, -n), y = n, colour = grupo)) +
  geom_point(size = 4.2) +
  geom_text(aes(label = n), vjust = -1.1, size = 3.8, colour = "grey20") +
  scale_colour_manual(values = c("Hepatocellular Carcinoma" = "#C0392B",
                                 "Hepatocholangiocarcinoma (Mixed)" = "#7D3C98",
                                 "Fibrolamellar Carcinoma" = "#F1C40F")) +
  scale_y_log10(labels = comma, breaks = c(1, 3, 10, 30, 100, 300)) +
  labs(title = "B. Subtipos histológicos (TCGA tumor)",
       subtitle = "Escala log10; n pequeno para mixed e fibrolamellar",
       x = NULL, y = "n (log10)", colour = "Histologia") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid.minor = element_blank())

p_grupos <- pA | pB
ggsave(file.path(dir_fig, "26_grupos_dados.png"), p_grupos,
       width = 11, height = 4.8, dpi = 300)

# -----------------------------------------------------------------------------
# 3. Marcadores por subtipo histológico (descritivo)
# -----------------------------------------------------------------------------
expr_file <- file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz")
expr_dt <- fread(expr_file, sep = "\t", data.table = FALSE)
rownames(expr_dt) <- expr_dt$gene
expr <- as.matrix(expr_dt[, setdiff(names(expr_dt), "gene")])
storage.mode(expr) <- "numeric"

marcadores <- intersect(c("GPC3", "AFP", "KRT19", "EPCAM", "BIRC5", "VIM"),
                        rownames(expr))
samples_hist <- tcga$sample[tcga$sample %in% colnames(expr)]
hist_long <- as.data.frame(t(expr[marcadores, samples_hist, drop = FALSE]))
hist_long$sample <- rownames(hist_long)
hist_long <- merge(hist_long, tcga[, c("sample", "histology")], by = "sample")
hist_long <- pivot_longer(hist_long, cols = all_of(marcadores),
                          names_to = "gene", values_to = "expressao")
hist_long$gene <- factor(hist_long$gene, levels = marcadores)

hist_summary <- hist_long %>%
  group_by(gene, histology) %>%
  summarise(
    n = n(),
    mediana = median(expressao),
    q1 = quantile(expressao, 0.25),
    q3 = quantile(expressao, 0.75),
    media = mean(expressao),
    .groups = "drop"
  )
write.csv(hist_summary, file.path(dir_tab, "22_histologia_marcadores.csv"),
          row.names = FALSE)

p_hist <- ggplot(hist_long, aes(x = histology, y = expressao, colour = histology)) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.7) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35,
               colour = "black", linewidth = 0.35, show.legend = FALSE) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("Hepatocellular Carcinoma" = "#C0392B",
                                 "Hepatocholangiocarcinoma (Mixed)" = "#7D3C98",
                                 "Fibrolamellar Carcinoma" = "#F1C40F")) +
  labs(title = "Marcadores por subtipo histológico — TCGA tumor",
       subtitle = "Cada ponto = uma amostra; barra preta = mediana. n mixed = 7, fibrolamellar = 3",
       x = NULL, y = "Expressão (log2 TPM + 0,001)", colour = "Histologia") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
ggsave(file.path(dir_fig, "27_histologia_marcadores.png"), p_hist,
       width = 10, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# 4. Resumo
# -----------------------------------------------------------------------------
resumo <- c(
  "GRUPOS ALTERNATIVOS — RESUMO",
  "============================",
  paste("Amostras totais:", nrow(meta_all), "| análise principal:", nrow(meta)),
  "",
  "Tipos de amostra:",
  paste(capture.output(print(tab_st, row.names = FALSE)), collapse = "\n"),
  "",
  "Histologia (TCGA tumor):",
  paste(capture.output(print(tab_hist, row.names = FALSE)), collapse = "\n"),
  "",
  "Sorologia viral (TCGA tumor):",
  paste(capture.output(print(tab_vir, row.names = FALSE)), collapse = "\n"),
  "",
  "Centros TCGA (TSS), top 10:",
  paste(capture.output(print(head(tab_tss, 10), row.names = FALSE)), collapse = "\n"),
  "",
  "Marcadores por histologia:",
  paste(capture.output(print(hist_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "Sessão R:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(resumo, file.path(dir_qc, "resumo_grupos_alternativos.txt"))
cat(paste(resumo, collapse = "\n"), "\n")
