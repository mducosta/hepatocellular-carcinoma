# =============================================================================
# Nanotec 2026 — GINAB/UENF
# Etapa 01 — Auditoria rigorosa dos dados (ANTES de qualquer análise)
# -----------------------------------------------------------------------------
# O QUE ESTE SCRIPT FAZ:
#   1. Lê o painel processado (41 genes x amostras de fígado).
#   2. Verifica: genes presentes vs ausentes; NA; variância zero; escala
#      (TPM linear vs log2); distribuição dos grupos; possível batch TCGA/GTEx.
#   3. Se os genes-chave (MMP9, MMP2, MMP14, TIMP1, COL1A1, COL4A1) estiverem
#      majoritariamente ausentes, GRAVA docs/AVISO_DADOS_INSUFICIENTES.txt e
#      PARA (não segue para análise).
#   4. Grava resumo em outputs/qc/auditoria.txt e docs/.
#
# Executar a partir da pasta nanotec2026/:
#   Rscript scripts/01_audit_data.R
# =============================================================================

cat("============================================================\n")
cat("NANOTEC 2026 - Etapa 01: Auditoria dos dados\n")
cat("============================================================\n")

base     <- "."
dir_raw  <- file.path(base, "dados", "raw")
dir_proc <- file.path(base, "dados", "processed")
dir_qc   <- file.path(base, "outputs", "qc")
dir_docs <- file.path(base, "docs")
dir.create(dir_qc, showWarnings = FALSE, recursive = TRUE)

panel <- c(
  "MMP1","MMP2","MMP3","MMP7","MMP8","MMP9","MMP10","MMP11","MMP12","MMP13",
  "MMP14","MMP15","MMP16",
  "TIMP1","TIMP2","TIMP3","TIMP4",
  "COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL6A1",
  "LAMC1","LAMA1","LAMB1",
  "FN1","SPP1","THBS1",
  "ITGB1","ITGAV","ITGA5","ITGA6",
  "PLAU","PLAUR","SERPINE1","VEGFA","TGFB1","CDH1","VIM"
)

# -----------------------------------------------------------------------------
# 1. Leitura do painel processado e da anotação
# -----------------------------------------------------------------------------
f_panel <- file.path(dir_proc, "liver_panel_41genes.tsv")
f_annot <- file.path(dir_proc, "samples_annotation.tsv")

if (!file.exists(f_panel)) {
  stop("Arquivo do painel não encontrado: ", f_panel,
       ". Rode scripts/00_download_data.R primeiro.")
}

cat("[LEITURA] painel:", f_panel, "\n")
mat <- read.delim(f_panel, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, quote = "")
rownames(mat) <- mat[[1]]          # 1ª coluna = símbolo gênico
mat <- mat[, -1, drop = FALSE]     # demais colunas = amostras

annot <- read.delim(f_annot, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, quote = "")

cat(sprintf("[LEITURA] matriz do painel: %d genes x %d amostras.\n",
            nrow(mat), ncol(mat)))

# -----------------------------------------------------------------------------
# 2. Genes presentes vs ausentes
# -----------------------------------------------------------------------------
genes_na_matriz <- rownames(mat)
presentes <- panel[panel %in% genes_na_matriz]
ausentes  <- panel[!panel %in% genes_na_matriz]

cat("\n[GENES] Presentes:", length(presentes), "/", length(panel), "\n")
if (length(ausentes) > 0) {
  cat("[GENES] AUSENTES:", paste(ausentes, collapse = ", "), "\n")
} else {
  cat("[GENES] Todos os genes do painel presentes.\n")
}

# -----------------------------------------------------------------------------
# 3. Verificação dos genes-chave (critério de parada)
# -----------------------------------------------------------------------------
genes_chave <- c("MMP9","MMP2","MMP14","TIMP1","COL1A1","COL4A1")
chave_presentes <- genes_chave[genes_chave %in% presentes]
cat("\n[GENES-CHAVE] Presentes:", length(chave_presentes), "/", length(genes_chave),
    "->", paste(chave_presentes, collapse = ", "), "\n")

if (length(chave_presentes) < 4) {
  aviso <- paste0(
    "AVISO DE DADOS INSUFICIENTES - Nanotec 2026\n",
    "============================================\n",
    "Data: ", Sys.time(), "\n\n",
    "Dos 6 genes-chave do protocolo (MMP9, MMP2, MMP14, TIMP1, COL1A1, COL4A1),\n",
    "apenas ", length(chave_presentes), " estao presentes no download:\n",
    "  Presentes: ", paste(chave_presentes, collapse = ", "), "\n",
    "  Ausentes : ", paste(setdiff(genes_chave, chave_presentes), collapse = ", "), "\n\n",
    "Motivo da parada: sem evidencia minima suficiente para uma analise de\n",
    "expressao diferencial robusta do eixo MMP9/MMPs/TIMPs na MEC do HCC.\n",
    "Nao e honesto prosseguir com menos da metade dos genes-chave.\n"
  )
  writeLines(aviso, file.path(dir_docs, "AVISO_DADOS_INSUFICIENTES.txt"))
  cat("\n[PARADA]", aviso)
  stop("Auditoria reprovada: genes-chave insuficientes. Veja docs/AVISO_DADOS_INSUFICIENTES.txt")
}
cat("[OK] Genes-chave suficientes. Auditoria prossegue.\n")

# -----------------------------------------------------------------------------
# 4. Valores faltantes e variância zero
# -----------------------------------------------------------------------------
mat_num <- as.matrix(mat)
storage.mode(mat_num) <- "double"

n_na <- sum(is.na(mat_num))
cat("\n[NA] Valores ausentes:", n_na, sprintf("(%.4f%%)", 100 * n_na / length(mat_num)), "\n")

var0 <- apply(mat_num, 1, var, na.rm = TRUE)
genes_var0 <- names(var0)[var0 == 0 | is.na(var0)]
cat("[VAR] Genes com variância zero:", length(genes_var0), "\n")
if (length(genes_var0) > 0) cat("       ->", paste(genes_var0, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 5. Escala dos dados (TPM linear vs log2)
# -----------------------------------------------------------------------------
cat("\n[ESCALA]\n")
cat("  min =", min(mat_num, na.rm = TRUE), "\n")
cat("  max =", max(mat_num, na.rm = TRUE), "\n")
cat("  mediana =", median(mat_num, na.rm = TRUE), "\n")
cat("  % valores > 100 =", round(100 * mean(mat_num > 100, na.rm = TRUE), 2), "%\n")
cat("  % valores < 1   =", round(100 * mean(mat_num < 1, na.rm = TRUE), 2), "%\n")

# Heurística de escala: TPM linear costuma ter max >> 20 e muitos valores > 1;
# log2(TPM+1) tem range tipicamente 0..~20.
parece_log2 <- max(mat_num, na.rm = TRUE) <= 25
cat("  Diagnóstico: dados", ifelse(parece_log2, "PODEM ser log2 (range <= 25).",
                                   "são TPM linear (max > 25)."), "\n")
cat("  -> Para análise usaremos log2(TPM+1) se forem TPM linear;\n",
    "     se já forem log2, usaremos como estão (verificado no script 02).\n")

# -----------------------------------------------------------------------------
# 6. Distribuição dos grupos
# -----------------------------------------------------------------------------
cat("\n[GRUPOS]\n")
if (all(annot[["sample"]] %in% colnames(mat))) {
  m <- match(colnames(mat), annot[["sample"]])
  grupos <- annot[["grupo"]][m]
} else {
  grupos <- rep(NA, ncol(mat))
  cat("  ATENÇÃO: nem todas as amostras da matriz estão na anotação!\n")
}
print(table(grupos, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 7. Possível batch effect TCGA vs GTEx (PCA)
# -----------------------------------------------------------------------------
cat("\n[BATCH] Gerando PCA do painel para inspeção de batch TCGA/GTEx...\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))

  # log2(TPM+1) se necessário
  X <- mat_num
  if (!parece_log2) X <- log2(X + 1)
  # genes sem variância (se houver) saem do PCA
  keep <- apply(X, 1, sd, na.rm = TRUE) > 0
  Xp <- t(X[keep, , drop = FALSE])
  pca <- prcomp(Xp, center = TRUE, scale. = TRUE)

  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], grupo = grupos)
  df[["estudo"]] <- ifelse(grepl("^GTEX", colnames(mat)), "GTEx", "TCGA")

  p1 <- ggplot(df, aes(PC1, PC2, color = grupo)) +
    geom_point(alpha = 0.7, size = 1.6) +
    labs(title = "PCA do painel ECM/MMP/TIMP — por grupo",
         x = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
         y = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2])) +
    theme_minimal(base_size = 12)
  ggsave(file.path(dir_qc, "PCA_painel_grupo.png"), p1, width = 7, height = 5, dpi = 150)

  p2 <- ggplot(df, aes(PC1, PC2, color = estudo)) +
    geom_point(alpha = 0.7, size = 1.6) +
    labs(title = "PCA do painel ECM/MMP/TIMP — por estudo (TCGA vs GTEx)",
         x = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
         y = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2])) +
    theme_minimal(base_size = 12)
  ggsave(file.path(dir_qc, "PCA_painel_estudo.png"), p2, width = 7, height = 5, dpi = 150)

  cat("[BATCH] PCA salvo em outputs/qc/ (PCA_painel_grupo.png e PCA_painel_estudo.png).\n")
} else {
  cat("[BATCH] ggplot2 indisponível; PCA pulado.\n")
}

# -----------------------------------------------------------------------------
# 8. Resumo da auditoria
# -----------------------------------------------------------------------------
linhas <- c(
  "RELATORIO DE AUDITORIA - Nanotec 2026 (Etapa 01)",
  paste0("Data: ", Sys.time()),
  "",
  paste0("Genes do painel (protocolo): ", length(panel)),
  paste0("Genes presentes na matriz: ", length(presentes)),
  paste0("Genes ausentes: ", ifelse(length(ausentes) > 0, paste(ausentes, collapse=", "), "nenhum")),
  paste0("Genes-chave presentes: ", paste(chave_presentes, collapse=", ")),
  "",
  paste0("Amostras na matriz: ", ncol(mat)),
  paste0("Grupos: ", paste(names(table(grupos)), table(grupos), sep="=", collapse=", ")),
  "",
  paste0("Valores ausentes (NA): ", n_na),
  paste0("Genes com variância zero: ", length(genes_var0)),
  paste0("Escala: ", ifelse(parece_log2, "possivelmente log2 (range <= 25)",
                            "TPM linear (max > 25)")),
  paste0("min=", min(mat_num, na.rm=TRUE), " max=", max(mat_num, na.rm=TRUE),
         " mediana=", median(mat_num, na.rm=TRUE)),
  "",
  "Decisão: auditoria APROVADA. Análise pode prosseguir (scripts 02-04)."
)
writeLines(linhas, file.path(dir_qc, "auditoria.txt"))
writeLines(linhas, file.path(dir_docs, "AUDITORIA_DADOS.txt"))

cat("\n[FIM] Auditoria concluída. Resumo salvo em outputs/qc/auditoria.txt e docs/AUDITORIA_DADOS.txt\n")
