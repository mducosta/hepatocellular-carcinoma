# =============================================================================
# etapa_04_volcano_plots.R
# -----------------------------------------------------------------------------
# Volcano plots para:
#   1. LIHC (TCGA): tumor primário vs tecido adjacente
#   2. GTEx: masculino vs feminino (comparação interna)
#   3. TCGA vs GTEx: tumor primário vs fígado normal
#   + comparações complementares (estágio, grau, sexo, idade, via LA, ComBat)
#
# Eixos:
#   x = logFC (log2 fold change)
#   y = -log10(FDR ajustado por BH)
#   Linhas: FDR = 0,05 e |logFC| = 1
# =============================================================================

set.seed(2026)  # reprodutibilidade dos rótulos/ggrepel

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_de  <- file.path(ROOT, "resultados", "tabelas", "expressao_diferencial")
dir_fig <- file.path(ROOT, "resultados", "figuras")
dir_tab <- file.path(ROOT, "resultados", "tabelas")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

FDR_CUT <- 0.05
LFC_CUT <- 1

cores <- c("Up" = "#C0392B", "Down" = "#2471A3", "NS" = "grey75")

fazer_volcano <- function(file_de, titulo, subtitulo, arquivo,
                          genes_rotular = NULL, n_rotulos = 10,
                          xlim = NULL, ylim = NULL) {
  d <- fread(file_de, sep = ",", data.table = FALSE)
  if (nrow(d) == 0) return(invisible(NULL))

  d$regulation <- factor(d$regulation, levels = c("Up", "Down", "NS"))
  d$neglog10 <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))

  # Rótulos: DEGs mais extremos (|logFC|) e, se solicitado, genes de interesse
  sig <- d[d$regulation != "NS", ]
  sig <- sig[order(-abs(sig$logFC)), ]
  lab <- head(sig, n_rotulos)
  if (!is.null(genes_rotular)) {
    extra <- d[d$gene %in% genes_rotular & !(d$gene %in% lab$gene), ]
    lab <- rbind(lab, extra)
  }
  if (nrow(lab) == 0) {
    lab <- head(d[order(d$adj.P.Val), ], n_rotulos)
  }

  n_up <- sum(d$regulation == "Up")
  n_down <- sum(d$regulation == "Down")
  sub <- paste0(subtitulo, " | Up = ", n_up, ", Down = ", n_down,
                " (FDR < 0,05 e |logFC| > 1)")

  p <- ggplot(d, aes(x = logFC, y = neglog10, colour = regulation)) +
    geom_point(size = 1.15, alpha = 0.65, stroke = 0) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
               linetype = "dashed", colour = "grey35", linewidth = 0.35) +
    geom_hline(yintercept = -log10(FDR_CUT),
               linetype = "dashed", colour = "grey35", linewidth = 0.35) +
    scale_colour_manual(values = cores, drop = FALSE,
                        labels = c("Up" = "Up", "Down" = "Down", "NS" = "Não significativo")) +
    labs(title = titulo, subtitle = sub,
         x = "log2 fold change (logFC)", y = "-log10(FDR ajustado)",
         colour = "Classe") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey30"),
      legend.position = "right"
    )

  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  p <- p + geom_text_repel(
    data = lab, aes(label = gene),
    size = 3, max.overlaps = 30, min.segment.length = 0,
    segment.colour = "grey50", segment.size = 0.25,
    box.padding = 0.35, point.padding = 0.2, seed = 2026
  )

  ggsave(file.path(dir_fig, paste0(arquivo, ".png")), p,
         width = 9, height = 7, dpi = 300)
  write.csv(lab[, c("gene", "logFC", "adj.P.Val", "regulation")],
            file.path(dir_tab, paste0("rotulos_", arquivo, ".csv")),
            row.names = FALSE)
  invisible(p)
}

# -----------------------------------------------------------------------------
# 1. Volcano principal — LIHC (TCGA): tumor vs adjacente
# -----------------------------------------------------------------------------
fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_TCGA_Adjacent_completo.csv"),
  "Volcano — LIHC (TCGA): tumor primário vs tecido adjacente",
  "Contraste: tumor − adjacente; limma com sexo e idade como covariáveis",
  "volcano_01_LIHC_TCGA_tumor_vs_adjacente",
  genes_rotular = c("AFP", "GPC3", "SPINK1", "MMP9", "COL1A1", "VIM", "CDH1", "BIRC5")
)

# -----------------------------------------------------------------------------
# 2. Volcano principal — GTEx: masculino vs feminino
# -----------------------------------------------------------------------------
fazer_volcano(
  file.path(dir_de, "expressao_diferencial_GTEx_Male_vs_Female_completo.csv"),
  "Volcano — GTEx (fígado normal): masculino vs feminino",
  "Comparação interna no GTEx; limma com idade como covariável",
  "volcano_02_GTEx_masculino_vs_feminino",
  genes_rotular = c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")
)

# -----------------------------------------------------------------------------
# 3. Volcano principal — TCGA vs GTEx: tumor vs fígado normal
# -----------------------------------------------------------------------------
fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_completo.csv"),
  "Volcano — TCGA-LIHC (tumor) vs GTEx (fígado normal)",
  "Contraste: tumor TCGA − fígado normal GTEx; limma com sexo e idade",
  "volcano_03_TCGA_tumor_vs_GTEx_normal",
  genes_rotular = c("AFP", "GPC3", "SPINK1", "MMP9", "COL1A1", "VIM", "CDH1",
                    "BIRC5", "EPCAM", "KRT19", "SPP1", "CXCL8")
)

# -----------------------------------------------------------------------------
# 4. Volcano complementares
# -----------------------------------------------------------------------------
fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_Adjacent_PAREADO_completo.csv"),
  "Volcano — LIHC (TCGA): tumor vs adjacente (análise pareada)",
  "Apenas pacientes com tumor e adjacente; paciente como fator de bloqueio",
  "volcano_04_LIHC_tumor_vs_adjacente_pareado",
  genes_rotular = c("AFP", "GPC3", "SPINK1", "MMP9", "COL1A1", "VIM")
)

fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Stage_III_IV_vs_I_II_completo.csv"),
  "Volcano — TCGA-LIHC: estágio III-IV vs I-II",
  "Apenas tumores primários; limma com sexo e idade",
  "volcano_05_TCGA_estagio_III_IV_vs_I_II"
)

fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Grade_G3_G4_vs_G1_G2_completo.csv"),
  "Volcano — TCGA-LIHC: grau G3-G4 vs G1-G2",
  "Apenas tumores primários; limma com sexo e idade",
  "volcano_06_TCGA_grau_G3_G4_vs_G1_G2"
)

fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Male_vs_Female_completo.csv"),
  "Volcano — TCGA-LIHC: masculino vs feminino",
  "Apenas tumores primários; limma com idade",
  "volcano_07_TCGA_masculino_vs_feminino"
)

fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Age65plus_vs_lt65_completo.csv"),
  "Volcano — TCGA-LIHC: idade ≥65 vs <65 anos",
  "Apenas tumores primários; limma com sexo",
  "volcano_08_TCGA_idade_65plus_vs_lt65"
)

fazer_volcano(
  file.path(dir_de, "expressao_diferencial_GTEx_Age60plus_vs_lt60_completo.csv"),
  "Volcano — GTEx: idade ≥60 vs <60 anos",
  "Fígado normal; limma com sexo",
  "volcano_09_GTEx_idade_60plus_vs_lt60"
)

fazer_volcano(
  file.path(dir_de, "expressao_diferencial_TCGA_Adjacent_vs_GTEx_Normal_completo.csv"),
  "Volcano — TCGA adjacente vs GTEx normal",
  "Comparação de tecidos normais de estudos diferentes (efeito de lote)",
  "volcano_10_TCGA_adjacente_vs_GTEx_normal"
)

fazer_volcano(
  file.path(dir_de, "expressao_diferencial_LApathway_TCGA_Tumor_vs_GTEx_Normal_completo.csv"),
  "Volcano — Via KEGG hsa05417 (Lipídios e Aterosclerose)",
  "Tumor TCGA vs fígado normal GTEx, restrito aos 213 genes da via",
  "volcano_11_LA_pathway_tumor_vs_GTEx"
)

# -----------------------------------------------------------------------------
# 4b. Volcano comparativo do ComBat (antes vs depois)
# -----------------------------------------------------------------------------
# O volcano "depois do ComBat" isolado é cientificamente quase vazio: apenas
# 1 gene permanece significativo. Para evitar uma figura visualmente vazia e
# mostrar o efeito da sobrecorreção, usamos dois painéis (antes/depois).
fazer_volcano_combat <- function(arquivo_antes, arquivo_depois, saida) {
  d1 <- fread(arquivo_antes, sep = ",", data.table = FALSE)
  d2 <- fread(arquivo_depois, sep = ",", data.table = FALSE)
  d1$etapa <- "Antes do ComBat"
  d2$etapa <- "Depois do ComBat (sem covariável)"
  d <- rbind(d1, d2)
  d$regulation <- factor(d$regulation, levels = c("Up", "Down", "NS"))
  d$etapa <- factor(d$etapa, levels = c("Antes do ComBat",
                                        "Depois do ComBat (sem covariável)"))
  d$neglog10 <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))

  # Rótulos: até 6 DEGs mais extremos por painel
  lab <- do.call(rbind, lapply(split(d, d$etapa), function(x) {
    s <- x[x$regulation != "NS", ]
    s <- s[order(-abs(s$logFC)), ]
    head(s, 6)
  }))

  n_antes <- sum(d1$regulation != "NS")
  n_depois <- sum(d2$regulation != "NS")

  p <- ggplot(d, aes(x = logFC, y = neglog10, colour = regulation)) +
    geom_point(size = 1.05, alpha = 0.6, stroke = 0) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
               linetype = "dashed", colour = "grey35", linewidth = 0.35) +
    geom_hline(yintercept = -log10(FDR_CUT),
               linetype = "dashed", colour = "grey35", linewidth = 0.35) +
    geom_text_repel(data = lab, aes(label = gene), size = 3,
                    max.overlaps = 30, min.segment.length = 0,
                    segment.colour = "grey50", segment.size = 0.25,
                    box.padding = 0.35, point.padding = 0.2, seed = 2026) +
    facet_wrap(~ etapa, scales = "free") +
    scale_colour_manual(values = cores, drop = FALSE,
                        labels = c("Up" = "Up", "Down" = "Down",
                                   "NS" = "Não significativo")) +
    labs(
      title = "Volcano — TCGA tumor vs GTEx normal: antes e depois do ComBat",
      subtitle = paste0("Antes: ", n_antes, " DEGs | Depois: ", n_depois,
                        " DEG — sobrecorreção esperada quando lote e condição são confundidos"),
      x = "log2 fold change (logFC)", y = "-log10(FDR ajustado)",
      colour = "Classe"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          strip.background = element_rect(fill = "grey92", colour = NA))

  ggsave(file.path(dir_fig, paste0(saida, ".png")), p,
         width = 12, height = 6, dpi = 300)
  write.csv(lab[, c("gene", "logFC", "adj.P.Val", "regulation", "etapa")],
            file.path(dir_tab, paste0("rotulos_", saida, ".csv")),
            row.names = FALSE)
  invisible(p)
}

fazer_volcano_combat(
  file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_completo.csv"),
  file.path(dir_de, "expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_ComBat_completo.csv"),
  "volcano_12_TCGA_vs_GTEx_ComBat"
)

cat("[VOLCANO] Figuras geradas em:", dir_fig, "\n")
# -----------------------------------------------------------------------------
# 5. Resumo dos top 5 genes por comparação (a partir dos arquivos TOP20)
# -----------------------------------------------------------------------------
top_files <- list.files(dir_de, pattern = "^top20_.*[.]csv$", full.names = TRUE)
if (length(top_files) > 0) {
  res_top <- do.call(rbind, lapply(top_files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    comp <- d$comparison[1]
    do.call(rbind, lapply(c("Up", "Down"), function(reg) {
      sub <- head(d[d$regulation == reg, ], 5)
      data.frame(
        comparacao = comp,
        regulacao = reg,
        top5_genes = paste(sub$gene, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
  }))
  write.csv(res_top, file.path(dir_de, "resumo_top5_genes.csv"), row.names = FALSE)
}

print(list.files(dir_fig, pattern = "^volcano", full.names = FALSE))
