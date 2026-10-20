# =============================================================================
# 04_ComBat_demonstracao.R
# -----------------------------------------------------------------------------
# Demonstração didática do ComBat (sva) no contexto TCGA-LIHC x GTEx.
#
# Cenário 1 (tentativa "correta"): ComBat com batch = estudo (TCGA/GTEx) e
#   mod = ~ condição (tumor/normal). Como a condição é completamente confundida
#   com o estudo (todo tumor é TCGA; todo normal é GTEx), o ComBat recusa:
#   "The covariate is confounded with batch!".
#
# Cenário 2 (demonstração de sobrecorreção): ComBat sem covariável biológica.
#   O lote é removido, mas também remove a diferença tumor-normal, pois ela é
#   indistinguível do lote. Mostramos a queda dos DEGs e o UMAP antes/depois.
#
# Este script NÃO substitui a análise principal; ele ilustra boas práticas e
# limitações. A análise principal (02) NÃO aplica ComBat por causa da
# confusão estrutural entre estudo e condição.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sva)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_proc <- file.path(ROOT, "dados", "processed")
dir_de   <- file.path(ROOT, "resultados", "tabelas", "DE")
dir_fig  <- file.path(ROOT, "resultados", "figuras")
dir_qc   <- file.path(ROOT, "resultados", "qc")
dir.create(dir_de, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

FDR_CUT <- 0.05
LFC_CUT <- 1

meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
              sep = "\t", data.table = FALSE)
meta$group <- factor(meta$group, levels = c("TCGA_Tumor", "TCGA_Adjacent", "GTEx_Normal"))
meta$sex <- factor(meta$sex, levels = c("female", "male"))

expr_dt <- fread(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"),
                 sep = "\t", data.table = FALSE)
rownames(expr_dt) <- expr_dt$gene
expr <- as.matrix(expr_dt[, setdiff(names(expr_dt), "gene")])
storage.mode(expr) <- "numeric"

meta_cb <- meta[meta$group %in% c("TCGA_Tumor", "GTEx_Normal"), ]
meta_cb$group_tg <- factor(ifelse(meta_cb$group == "TCGA_Tumor", "TCGA_Tumor", "GTEx_Normal"),
                           levels = c("GTEx_Normal", "TCGA_Tumor"))
meta_cb <- meta_cb[complete.cases(meta_cb[, c("sex", "age")]), ]
X <- expr[, meta_cb$sample, drop = FALSE]
batch <- factor(meta_cb$study, levels = c("TCGA", "GTEX"))

log_lines <- c("DEMONSTRAÇÃO COMBAT", "===================",
               paste("Amostras:", ncol(X), "| Genes:", nrow(X)),
               paste("TCGA tumor:", sum(meta_cb$group_tg == "TCGA_Tumor"),
                     "| GTEx normal:", sum(meta_cb$group_tg == "GTEx_Normal")))

# -----------------------------------------------------------------------------
# 1. ComBat com preservação da condição (deve falhar por confusão)
# -----------------------------------------------------------------------------
mod_bio <- model.matrix(~ group_tg, data = meta_cb)
erro_confusao <- tryCatch({
  invisible(sva::ComBat(dat = X, batch = batch, mod = mod_bio, par.prior = TRUE))
  NA_character_
}, error = function(e) e$message)

log_lines <- c(log_lines, "",
               "Cenário 1 — ComBat com mod = ~ grupo:",
               paste("Resultado:", ifelse(is.na(erro_confusao),
                                           "executado (inesperado)",
                                           erro_confusao)))

# -----------------------------------------------------------------------------
# 2. ComBat sem covariável (sobrecorreção)
# -----------------------------------------------------------------------------
set.seed(2026)
X_cor <- sva::ComBat(dat = X, batch = batch, mod = NULL, par.prior = TRUE)
log_lines <- c(log_lines, "",
               "Cenário 2 — ComBat sem covariável biológica:",
               "Lote removido, mas a condição tumor/normal é confundida com o lote.")

# Função limma
rodar_limma <- function(Xmat, meta_df, nome) {
  design <- model.matrix(~ 0 + group_tg + sex + age, data = meta_df)
  colnames(design) <- sub("group_tg", "group", colnames(design), fixed = TRUE)
  fit <- lmFit(Xmat, design)
  cvec <- rep(0, ncol(design)); names(cvec) <- colnames(design)
  cvec["groupGTEx_Normal"] <- -1
  cvec["groupTCGA_Tumor"] <-  1
  fit2 <- contrasts.fit(fit, cvec)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt$comparison <- nome
  tt$direction <- ifelse(tt$logFC > 0, "TCGA_Tumor > GTEx_Normal", "GTEx_Normal > TCGA_Tumor")
  tt$regulation <- ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC >  LFC_CUT, "Up",
                   ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC < -LFC_CUT, "Down", "NS"))
  tt[, c("gene", "comparison", "direction", "logFC", "AveExpr",
         "t", "P.Value", "adj.P.Val", "B", "regulation")]
}

res_antes <- rodar_limma(X, meta_cb, "TCGA_Tumor_vs_GTEx_Normal")
res_depois <- rodar_limma(X_cor, meta_cb, "TCGA_Tumor_vs_GTEx_Normal_ComBat_semCovariavel")

write.csv(res_antes, file.path(dir_de, "DEG_TCGA_Tumor_vs_GTEx_Normal_completo.csv"),
          row.names = FALSE)
write.csv(res_depois, file.path(dir_de, "DEG_TCGA_Tumor_vs_GTEx_Normal_ComBat_completo.csv"),
          row.names = FALSE)
# Mantém também o nome usado pelo script de volcano
write.csv(res_depois, file.path(dir_de, "DEG_TCGA_Tumor_vs_GTEx_Normal_ComBat_semCovariavel.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Comparação antes/depois
# -----------------------------------------------------------------------------
comp <- data.frame(
  metrica = c("DEGs antes (FDR<0.05, |logFC|>1)",
              "DEGs depois do ComBat sem covariável",
              "Up antes", "Up depois", "Down antes", "Down depois",
              "DEGs em comum", "Jaccard dos DEGs",
              "Correlação de logFC (Spearman)"),
  valor = c(
    sum(res_antes$regulation != "NS"),
    sum(res_depois$regulation != "NS"),
    sum(res_antes$regulation == "Up"),
    sum(res_depois$regulation == "Up"),
    sum(res_antes$regulation == "Down"),
    sum(res_depois$regulation == "Down"),
    length(intersect(res_antes$gene[res_antes$regulation != "NS"],
                     res_depois$gene[res_depois$regulation != "NS"])),
    round(length(intersect(res_antes$gene[res_antes$regulation != "NS"],
                           res_depois$gene[res_depois$regulation != "NS"])) /
            max(1, length(union(res_antes$gene[res_antes$regulation != "NS"],
                                res_depois$gene[res_depois$regulation != "NS"]))), 3),
    {
      common_genes <- intersect(res_antes$gene, res_depois$gene)
      round(cor(res_antes$logFC[match(common_genes, res_antes$gene)],
                res_depois$logFC[match(common_genes, res_depois$gene)],
                method = "spearman", use = "complete.obs"), 3)
    }
  )
)
write.csv(comp, file.path(dir_de, "COMPARACAO_ComBat_antes_depois.csv"), row.names = FALSE)
log_lines <- c(log_lines, "", "Comparação antes/depois:",
               paste(capture.output(print(comp, row.names = FALSE)), collapse = "\n"))

# -----------------------------------------------------------------------------
# 4. UMAP antes/depois
# -----------------------------------------------------------------------------
vars <- apply(X, 1, var, na.rm = TRUE)
topv <- names(sort(vars, decreasing = TRUE))[seq_len(min(2000, length(vars)))]
set.seed(2026)
umap_b <- uwot::umap(t(X[topv, ]), n_neighbors = 30, min_dist = 0.25,
                     metric = "cosine", scale = TRUE, pca = 50,
                     seed = 2026, n_threads = 1, n_sgd_threads = 1,
                     verbose = FALSE)
umap_a <- uwot::umap(t(X_cor[topv, ]), n_neighbors = 30, min_dist = 0.25,
                     metric = "cosine", scale = TRUE, pca = 50,
                     seed = 2026, n_threads = 1, n_sgd_threads = 1,
                     verbose = FALSE)
umap_df <- rbind(
  data.frame(sample = colnames(X), UMAP1 = umap_b[, 1], UMAP2 = umap_b[, 2],
             etapa = "Antes do ComBat"),
  data.frame(sample = colnames(X_cor), UMAP1 = umap_a[, 1], UMAP2 = umap_a[, 2],
             etapa = "Depois do ComBat (sem covariável)")
)
umap_df <- merge(umap_df, meta_cb[, c("sample", "group_tg", "study")], by = "sample")
write.csv(umap_df, file.path(dir_qc, "UMAP_ComBat_coordenadas.csv"), row.names = FALSE)

p <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = group_tg, shape = study)) +
  geom_point(size = 1.7, alpha = 0.75) +
  facet_wrap(~ etapa, scales = "free") +
  scale_colour_manual(values = c("TCGA_Tumor" = "#C0392B",
                                 "GTEx_Normal" = "#2471A3")) +
  labs(title = "ComBat: antes e depois da correção de lote TCGA x GTEx",
       subtitle = "Demonstração de sobrecorreção quando lote e condição são confundidos",
       x = "UMAP 1", y = "UMAP 2", colour = "Grupo", shape = "Estudo") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(dir_fig, "13_UMAP_ComBat_antes_depois.png"), p,
       width = 10, height = 5, dpi = 300)

# Violino + boxplot de um marcador (AFP, se presente) antes/depois
if ("AFP" %in% rownames(X)) {
  box_df <- data.frame(
    sample = colnames(X),
    grupo = meta_cb$group_tg,
    Antes = as.numeric(X["AFP", ]),
    Depois = as.numeric(X_cor["AFP", ])
  )
  box_long <- rbind(
    data.frame(sample = box_df$sample, grupo = box_df$grupo,
               etapa = "Antes do ComBat", expressao = box_df$Antes),
    data.frame(sample = box_df$sample, grupo = box_df$grupo,
               etapa = "Depois do ComBat (sem covariável)", expressao = box_df$Depois)
  )
  p2 <- ggplot(box_long, aes(x = grupo, y = expressao, fill = grupo)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.55, colour = NA) +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9,
                 colour = "grey20") +
    facet_wrap(~ etapa) +
    scale_fill_manual(values = c("TCGA_Tumor" = "#C0392B",
                                 "GTEx_Normal" = "#2471A3")) +
    labs(title = "Exemplo: expressão de AFP antes e depois do ComBat (violino + boxplot)",
         subtitle = "A remoção do lote também remove a diferença biológica esperada",
         x = NULL, y = "log2(TPM + 0,001)", fill = "Grupo") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(dir_fig, "14_AFP_ComBat_antes_depois.png"), p2,
         width = 8, height = 4.5, dpi = 300)
}

writeLines(log_lines, file.path(dir_qc, "resumo_ComBat_demonstracao.txt"))
cat(paste(log_lines, collapse = "\n"), "\n")
