# ===================================================================
# ANÁLISE DE EXPRESSÃO DIFERENCIAL — LIHC vs Normal (GTEx)
# Via Lipid and Atherosclerosis (KEGG hsa05417)
# ===================================================================
# PORTÁTIL · REPRODUTÍVEL · AUTODOCUMENTADO
# Executar: Rscript pipeline_hepato.R
# ===================================================================
# DECLARAÇÃO DE USO DE INTELIGÊNCIA ARTIFICIAL (Portaria CNPq nº 2.664/2026):
# "Este trabalho contou com o uso de ferramentas de inteligência artificial
#  generativa (ChatGPT-5.5, OpenAI; DeepSeek-V4-Pro, DeepSeek) como suporte
#  para processamento de dados transcriptômicos, revisão de código em
#  linguagem R e editoração científica, sob supervisão e validação integral
#  dos autores, que assumem total responsabilidade pelo conteúdo final."
# ===================================================================

# ------------------------------------------------------------------
# 0) CONFIGURAÇÃO INICIAL — PORTABILIDADE
# ------------------------------------------------------------------

# Suprimir Rplots.pdf
pdf(NULL)

# 0.1 Detectar raiz do projeto automaticamente
#     (procura pelo diretório que contém este script ou o arquivo de dados)
find_project_root <- function() {
  candidates <- c(getwd())
  # Tentar obter diretório do script (funciona com Rscript e source)
  tryCatch({
    script_path <- sys.frame(1)$ofile
    if (is.character(script_path)) candidates <- c(candidates, dirname(script_path))
  }, error = function(e) {})
  # Também tentar commandArgs
  tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- sub("--file=", "", file_arg[1])
      candidates <- c(candidates, dirname(script_path))
    }
  }, error = function(e) {})
  candidates <- c(candidates, ".")
  candidates <- unique(candidates)
  for (d in candidates) {
    d <- normalizePath(d, mustWork = FALSE)
    data_patterns <- c("liver\\.(tsv|csv|xlsx)$",
                       "liver_lip_aterosclerose\\.(tsv|csv|xlsx)$",
                       "liver_lip_aterosclero\\.(tsv|csv|xlsx)$")
    # Buscar também na subpasta data/ (estrutura padrão do repositório)
    search_dirs <- unique(c(d, file.path(d, "data")))
    for (sd in search_dirs) {
      for (pat in data_patterns) {
        if (length(list.files(sd, pattern = pat, ignore.case = TRUE)) > 0) {
          return(d)
        }
      }
    }
  }
  return(getwd())
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
cat(sprintf("PROJECT_ROOT: %s\n", PROJECT_ROOT))

# 0.2 Criar estrutura de diretórios de saída
dirs_out <- c("results", "results/qc", "results/deg", "results/volcano",
              "results/ppi", "results/enrichment", "results/audit",
              "results/figures", "results/tables")
for (d in dirs_out) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# 0.3 Log pipeline
LOG_FILE <- file.path(PROJECT_ROOT, "results/audit/pipeline_log.csv")
log_entry <- function(stage, status, detail) {
  line <- sprintf('"%s","%s","%s","%s"',
                  format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"),
                  stage, status, detail)
  cat(line, "\n", sep = "")
  write(line, file = LOG_FILE, append = TRUE)
}
write('"timestamp","stage","status","detail"', file = LOG_FILE)
log_entry("START", "OK", paste("R version:", R.version.string))
log_entry("PROJECT_ROOT", "OK", PROJECT_ROOT)

# 0.4 Timestamp global para benchmark
T0 <- Sys.time()

# ------------------------------------------------------------------
# 1) PACOTES
# ------------------------------------------------------------------
required_packages <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel",
  "limma", "edgeR", "igraph", "ggraph", "pheatmap", "scales",
  "clusterProfiler", "org.Hs.eg.db", "KEGGREST", "STRINGdb",
  "readr", "readxl", "rio", "uwot", "fgsea", "msigdbr", "GSVA"
)

log_entry("PACKAGES", "INFO", "Instalando/carregando pacotes...")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(edgeR)
  library(igraph)
  library(ggraph)
  library(pheatmap)
  library(scales)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(KEGGREST)
  library(STRINGdb)
  library(readr)
  library(readxl)
  library(rio)
  library(uwot)
  library(fgsea)
  library(msigdbr)
  library(GSVA)
})

log_entry("PACKAGES", "OK", paste(length(required_packages), "pacotes carregados"))

# ------------------------------------------------------------------
# 2) LOCALIZAR ARQUIVO DE DADOS (PORTÁTIL)
# ------------------------------------------------------------------
data_patterns <- c(
  "liver.tsv",
  "liver.csv",
  "liver.xlsx",
  "liver_lip_aterosclerose.tsv",
  "liver_lip_aterosclerose.csv",
  "liver_lip_aterosclerose.xlsx",
  "liver_lip_aterosclero.tsv",
  "liver_lip_aterosclero.csv",
  "liver_lip_aterosclero.xlsx"
)

# Buscar o arquivo de dados na raiz e na subpasta data/
data_search_dirs <- unique(c(PROJECT_ROOT, file.path(PROJECT_ROOT, "data")))
data_files_found <- unlist(lapply(data_search_dirs, function(d) {
  list.files(d, pattern = "\\.(tsv|csv|xlsx)$", ignore.case = TRUE,
             full.names = TRUE)
}))
data_files_found <- data_files_found[
  tolower(basename(data_files_found)) %in% tolower(data_patterns) |
  grepl("^liver", tolower(basename(data_files_found)))
]

if (length(data_files_found) == 0) {
  log_entry("DATA_FILE", "FAIL", "Nenhum arquivo de dados encontrado")
  stop("ARQUIVO DE DADOS NÃO ENCONTRADO. Coloque 'liver_lip_aterosclerose.tsv' na pasta do projeto.")
}

DATA_FILE <- data_files_found[1]
log_entry("DATA_FILE", "OK", DATA_FILE)

# 2.1 Importar de acordo com a extensão
ext <- tolower(tools::file_ext(DATA_FILE))
log_entry("IMPORT", "INFO", paste("Extensão:", ext))

if (ext == "tsv") {
  liver_raw <- readr::read_tsv(DATA_FILE, show_col_types = FALSE, name_repair = "minimal")
} else if (ext == "csv") {
  liver_raw <- readr::read_csv(DATA_FILE, show_col_types = FALSE, name_repair = "minimal")
} else if (ext %in% c("xls", "xlsx")) {
  liver_raw <- readxl::read_excel(DATA_FILE, .name_repair = "minimal")
} else {
  stop("Formato não suportado: ", ext)
}

liver <- as.data.frame(liver_raw)
log_entry("IMPORT", "OK", paste("Dimensão bruta:", paste(dim(liver), collapse = " x ")))

# ------------------------------------------------------------------
# 3) NORMALIZAR NOMES DE COLUNAS
# ------------------------------------------------------------------
# Remove underscores iniciais (Xena usa _primary_site, _sample_type, etc.)
names(liver) <- sub("^_+", "", names(liver))

# Remove coluna duplicada 'samples' se existir (redundante com 'sample')
if ("samples" %in% colnames(liver)) liver$samples <- NULL

# ------------------------------------------------------------------
# 4) IDENTIFICAR METADADOS vs GENES (EXPRESSÃO)
# ------------------------------------------------------------------
# As colunas de metadados esperadas da Xena (após limpeza de underscore)
meta_cols_expected <- c("sample", "primary_site", "sample_type", "study",
                         "TCGA_GTEX_main_category", "detailed_category",
                         "primary disease or tissue")

# Identificar quais colunas são realmente metadados (presentes no arquivo)
meta_cols_detected <- base::intersect(meta_cols_expected, colnames(liver))
gene_cols <- setdiff(colnames(liver), meta_cols_detected)

# Também remover colunas clínicas adicionais que podem estar presentes
clinical_keywords <- c("^OS$", "^OS\\.time$", "^DSS$", "^DSS\\.time$",
                       "^DFI$", "^DFI\\.time$", "^PFI$", "^PFI\\.time$",
                       "^_?gender$", "^RFS$", "^RFS\\.time$", "^DSS$",
                       "^age", "^stage", "^grade", "^vital")
for (kw in clinical_keywords) {
  extra_clinical <- grep(kw, gene_cols, value = TRUE, ignore.case = TRUE)
  if (length(extra_clinical) > 0) {
    meta_cols_detected <- c(meta_cols_detected, extra_clinical)
    gene_cols <- setdiff(gene_cols, extra_clinical)
  }
}

log_entry("COLUMNS", "OK", paste("meta:", length(meta_cols_detected),
                                  "genes:", length(gene_cols)))

# ------------------------------------------------------------------
# 5) VERIFICAÇÃO CRÍTICA: DADOS DE EXPRESSÃO GÊNICA PRESENTES?
# ------------------------------------------------------------------
if (length(gene_cols) < 10) {
  log_entry("EXPR_DATA", "FAIL",
            paste("Apenas", length(gene_cols),
                  "colunas candidatas a genes. Dados de expressão gênica AUSENTES."))
  cat("\n")
  cat("============================================================\n")
  cat(" ERRO CRÍTICO: ARQUIVO NÃO CONTÉM DADOS DE EXPRESSÃO GÊNICA\n")
  cat("============================================================\n")
  cat("\n")
  cat("O arquivo '", basename(DATA_FILE), "' contém apenas metadados clínicos.\n", sep = "")
  cat(length(meta_cols_detected), " colunas de metadados e ", length(gene_cols),
      " colunas adicionais, mas nenhuma coluna de expressão gênica.\n\n", sep = "")
  cat("PARA CORRIGIR:\n")
  cat("1. Acesse https://xenabrowser.net/\n")
  cat("2. Adicione os datasets:\n")
  cat("   - TCGA Liver Hepatocellular Carcinoma (LIHC) - gene expression RNAseq\n")
  cat("   - GTEx Liver - gene expression RNAseq\n")
  cat("3. Selecione 'Gene Expression' como visualização\n")
  cat("4. Baixe o arquivo TSV completo (File > Download)\n")
  cat("5. Substitua o arquivo atual pelo novo download\n\n")
  cat("Colunas de metadados encontradas:\n")
  cat(paste("  -", meta_cols_detected), sep = "\n")
  cat("\n")
  cat("Gerando relatório de auditoria mesmo sem dados de expressão...\n")
  
  # Salvar sessionInfo antes de parar
  writeLines(capture.output(sessionInfo()),
             file.path(PROJECT_ROOT, "results/audit/sessionInfo.txt"))
  
  # Benchmark mesmo com falha
  T1 <- Sys.time()
  bench <- data.frame(
    metric = "tempo_total_execucao_seg",
    value  = round(as.numeric(difftime(T1, T0, units = "secs")), 2),
    status = "FALHA_PARCIAL"
  )
  rio::export(bench, file.path(PROJECT_ROOT, "results/audit/benchmark_pipeline.csv"))
  
  log_entry("FATAL", "FAIL", "Dados de expressão gênica ausentes. Pipeline abortado.")
  
  # Gerar relatório de erro
  sink(file.path(PROJECT_ROOT, "results/audit/AUDITORIA_REANALISE_LIHC.md"))
  cat("# AUDITORIA DE REANÁLISE LIHC\n\n")
  cat("## STATUS: NÃO EXECUTADO — DADOS INCOMPLETOS\n\n")
  cat("**Data:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("### Problema identificado\n\n")
  cat("O arquivo de dados não contém colunas de expressão gênica.\n")
  cat("Possui apenas ", length(meta_cols_detected), " colunas de metadados.\n\n")
  cat("### Ação necessária\n\n")
  cat("Fazer download do dataset completo do UCSC Xena incluindo expressão gênica.\n")
  sink()
  
  stop("DADOS DE EXPRESSÃO GÊNICA AUSENTES. Consulte o relatório em results/audit/.")
}

# ------------------------------------------------------------------
# 6) FILTRAR APENAS AMOSTRAS DE FÍGADO (LIHC + GTEx Liver)
# ------------------------------------------------------------------
# Detectar a coluna que identifica o tipo de tecido
category_col <- "TCGA_GTEX_main_category"
if (!(category_col %in% colnames(liver))) {
  # tentar outras colunas: 'primary disease or tissue', 'detailed_category'
  if ("detailed_category" %in% colnames(liver)) {
    category_col <- "detailed_category"
  } else if ("primary disease or tissue" %in% colnames(liver)) {
    category_col <- "primary disease or tissue"
  }
}

liver_categories <- unique(liver[[category_col]])
cat("Categorias encontradas na coluna '", category_col, "':\n", sep = "")
print(head(liver_categories, 30))

# Filtrar para Liver/LIHC
liver_mask <- grepl("Liver|LIHC|Hepatocellular", liver[[category_col]], ignore.case = TRUE)
liver_sub <- liver[liver_mask, , drop = FALSE]
log_entry("FILTER_LIVER", "OK", paste("Amostras fígado:", nrow(liver_sub)))

if (nrow(liver_sub) < 50) {
  log_entry("FILTER_LIVER", "FAIL", "Menos de 50 amostras hepáticas. Verifique o arquivo.")
  stop("Poucas amostras hepáticas detectadas.")
}

# ------------------------------------------------------------------
# 7) SEPARAR METADADOS E EXPRESSÃO
# ------------------------------------------------------------------
meta <- liver_sub[, meta_cols_detected, drop = FALSE]
expr_cols <- base::intersect(gene_cols, colnames(liver_sub))

# Verificar se as colunas de gene são numéricas
gene_is_numeric <- sapply(liver_sub[expr_cols], is.numeric)
non_numeric_genes <- names(gene_is_numeric)[!gene_is_numeric]

if (length(non_numeric_genes) > 0) {
  cat("Aviso:", length(non_numeric_genes), "colunas de 'gene' não são numéricas. Convertendo...\n")
  for (col in non_numeric_genes) {
    liver_sub[[col]] <- as.numeric(as.character(liver_sub[[col]]))
  }
}

expr <- as.matrix(liver_sub[, expr_cols, drop = FALSE])
storage.mode(expr) <- "numeric"

cat(sprintf("Matriz de expressão: %d amostras x %d genes\n", nrow(expr), ncol(expr)))
log_entry("MATRIX", "OK", paste(nrow(expr), "x", ncol(expr)))
log_entry("PROGRESS", "INFO", "Matriz de expressão montada. Iniciando definição de condições...")

# ------------------------------------------------------------------
# 8) DEFINIR CONDIÇÕES (Normal vs LIHC) — COM FILTRO POR SAMPLE_TYPE
# ------------------------------------------------------------------
log_entry("PROGRESS", "INFO", "Definindo condições com filtro sample_type...")

# Verificar distribuição de sample_type por categoria
cat("\nDistribuição detalhada TCGA_GTEX_main_category × sample_type:\n")
print(table(meta$TCGA_GTEX_main_category, meta$sample_type))

# Definir condição com filtro explícito de sample_type (P1 da auditoria)
if ("TCGA_GTEX_main_category" %in% colnames(meta)) {
  meta$condition <- dplyr::case_when(
    grepl("GTEX.*Liver", meta$TCGA_GTEX_main_category, ignore.case = TRUE) &
      grepl("Normal", meta$sample_type, ignore.case = TRUE) ~ "Normal",
    grepl("TCGA.*Liver.*Hepatocellular|LIHC", meta$TCGA_GTEX_main_category,
          ignore.case = TRUE) &
      grepl("Primary", meta$sample_type, ignore.case = TRUE) ~ "LIHC",
    TRUE ~ NA_character_
  )
} else if ("study" %in% colnames(meta)) {
  meta$condition <- dplyr::case_when(
    meta$study == "GTEX" & grepl("Normal", meta$sample_type, ignore.case = TRUE) ~ "Normal",
    meta$study == "TCGA" & grepl("Primary", meta$sample_type, ignore.case = TRUE) ~ "LIHC",
    TRUE ~ NA_character_
  )
} else {
  stop("Não foi possível identificar condições Normal/LIHC.")
}

n_excluded <- sum(is.na(meta$condition))
cat(sprintf("Amostras excluídas por sample_type inadequado: %d\n", n_excluded))
log_entry("SAMPLE_TYPE_FILTER", "OK", paste(n_excluded, "amostras excluídas por sample_type"))

meta$condition <- factor(meta$condition, levels = c("Normal", "LIHC"))

# Remover amostras sem condição definida
keep <- !is.na(meta$condition)
meta <- meta[keep, , drop = FALSE]
expr <- expr[keep, , drop = FALSE]

stopifnot(sum(meta$condition == "Normal", na.rm = TRUE) >= 2)
stopifnot(sum(meta$condition == "LIHC", na.rm = TRUE) >= 2)

# Identificar batch/coorte
if ("study" %in% colnames(meta)) {
  cat("Distribuição por estudo:\n")
  print(table(meta$condition, meta$study))
} else {
  cat("Distribuição por condição:\n")
  print(table(meta$condition))
}

log_entry("CONDITIONS", "OK", paste("Normal:", sum(meta$condition == "Normal"),
                                     "LIHC:", sum(meta$condition == "LIHC")))

# ------------------------------------------------------------------
# 9) CONTROLE DE QUALIDADE PRÉ-ANÁLISE
# ------------------------------------------------------------------
log_entry("QC_START", "OK", "Iniciando controle de qualidade pré-análise")

qc_report <- list()

# 9.1 Valores ausentes
na_count <- sum(is.na(expr))
na_pct <- mean(is.na(expr)) * 100
qc_report$na_count <- na_count
qc_report$na_pct <- na_pct
log_entry("QC_NA", if (na_count == 0) "OK" else "WARN",
          paste("NA count:", na_count, sprintf("(%.2f%%)", na_pct)))

# 9.2 Escala dos dados
expr_range <- range(expr, na.rm = TRUE)
is_integer_like <- all(abs(expr - round(expr)) < 1e-10, na.rm = TRUE)
is_log_scale <- expr_range[1] >= 0 && expr_range[2] < 25
qc_report$expr_range <- expr_range
qc_report$is_integer_like <- is_integer_like
qc_report$is_log_scale <- is_log_scale
log_entry("QC_SCALE", "INFO",
          sprintf("Range: %.2f - %.2f | Integer-like: %s | Log-scale: %s",
                  expr_range[1], expr_range[2], is_integer_like, is_log_scale))

# 9.3 Genes duplicados
dup_genes <- names(which(table(colnames(expr)) > 1))
qc_report$dup_genes <- length(dup_genes)
if (length(dup_genes) > 0) {
  log_entry("QC_DUP_GENES", "WARN", paste(length(dup_genes), "genes duplicados"))
  # Resolver: média das duplicatas
  for (g in dup_genes) {
    idx <- which(colnames(expr) == g)
    expr[, idx[1]] <- rowMeans(expr[, idx, drop = FALSE], na.rm = TRUE)
    expr <- expr[, -idx[-1], drop = FALSE]
  }
} else {
  log_entry("QC_DUP_GENES", "OK", "Sem genes duplicados")
}

# 9.4 Amostras duplicadas
dup_samples <- names(which(table(meta$sample) > 1))
qc_report$dup_samples <- length(dup_samples)
log_entry("QC_DUP_SAMPLES", if (length(dup_samples) == 0) "OK" else "WARN",
          paste(length(dup_samples), "amostras duplicadas"))

# 9.5 Genes com variância zero (remover)
gene_var <- apply(expr, 2, var, na.rm = TRUE)
zero_var_genes <- colnames(expr)[gene_var == 0 | is.na(gene_var)]
qc_report$zero_var_genes <- length(zero_var_genes)
if (length(zero_var_genes) > 0) {
  expr <- expr[, !colnames(expr) %in% zero_var_genes, drop = FALSE]
  log_entry("QC_ZERO_VAR", "WARN", paste(length(zero_var_genes), "genes com variância zero removidos"))
} else {
  log_entry("QC_ZERO_VAR", "OK", "Sem genes de variância zero")
}

# 9.6 Verificação de batch/coorte
log_entry("QC_BATCH", "INFO", "Verificando efeito batch TCGA vs GTEx")

# 9.7 PCA pré-análise
cat("\n--- PCA Pré-Análise ---\n")
expr_t <- t(expr)  # genes x amostras para PCA
expr_t_scaled <- expr_t[apply(expr_t, 1, var, na.rm = TRUE) > 0, , drop = FALSE]

if (nrow(expr_t_scaled) >= 2 && ncol(expr_t_scaled) >= 3) {
  pca_res <- prcomp(t(expr_t_scaled), scale. = TRUE)
  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Condition = meta$condition
  )
  
  if ("study" %in% colnames(meta)) pca_df$Study <- meta$study
  
  pc1_var <- round(summary(pca_res)$importance[2, 1] * 100, 1)
  pc2_var <- round(summary(pca_res)$importance[2, 2] * 100, 1)
  
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95, linewidth = 0.8) +
    labs(
      title = "PCA Pré-Análise — Todos os Genes",
      subtitle = paste(nrow(meta), "amostras |", nrow(expr_t_scaled), "genes"),
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    ) +
    scale_color_manual(values = c(Normal = "#2E86AB", LIHC = "#A23B72")) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  # Se houver batch/study, gerar também PCA colorido por estudo
  if ("study" %in% colnames(meta)) {
    p_pca_study <- ggplot(pca_df, aes(PC1, PC2, color = Study, shape = Condition)) +
      geom_point(size = 3, alpha = 0.7) +
      stat_ellipse(aes(group = interaction(Condition, Study)), level = 0.95, linewidth = 0.5) +
      labs(
        title = "PCA Pré-Análise — Colorido por Estudo (Batch)",
        subtitle = "Verificação de separação TCGA vs GTEx",
        x = paste0("PC1 (", pc1_var, "%)"),
        y = paste0("PC2 (", pc2_var, "%)")
      ) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(PROJECT_ROOT, "results/qc/PCA_pre_analysis_batch.png"),
           p_pca_study, width = 9, height = 7, dpi = 300)
  }
  
  ggsave(file.path(PROJECT_ROOT, "results/qc/PCA_pre_analysis.png"),
         p_pca, width = 8, height = 6, dpi = 300)
  log_entry("QC_PCA", "OK", paste("PC1:", pc1_var, "% PC2:", pc2_var, "%"))
} else {
  log_entry("QC_PCA", "SKIP", "Dados insuficientes para PCA")
}

# 9.8 UMAP pré-análise
cat("\n--- UMAP Pré-Análise ---\n")
if (nrow(expr_t_scaled) >= 2 && ncol(expr_t_scaled) >= 5) {
  set.seed(4721)
  umap_res <- uwot::umap(t(expr_t_scaled), n_neighbors = min(15, nrow(meta) - 1),
                          min_dist = 0.3, metric = "euclidean")
  umap_df <- data.frame(
    UMAP1 = umap_res[, 1],
    UMAP2 = umap_res[, 2],
    Condition = meta$condition
  )
  
  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Condition)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = "UMAP Pré-Análise — Todos os Genes",
      subtitle = paste(nrow(meta), "amostras |", nrow(expr_t_scaled), "genes"),
      x = "UMAP1", y = "UMAP2"
    ) +
    scale_color_manual(values = c(Normal = "#2E86AB", LIHC = "#A23B72")) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(PROJECT_ROOT, "results/qc/UMAP_pre_analysis.png"),
         p_umap, width = 8, height = 6, dpi = 300)
  log_entry("QC_UMAP", "OK", "UMAP gerado com sucesso")
} else {
  log_entry("QC_UMAP", "SKIP", "Dados insuficientes para UMAP")
}

# 9.9 Distribuição de amostras
sample_dist <- as.data.frame(table(meta$condition))
if ("study" %in% colnames(meta)) {
  sample_dist_study <- as.data.frame(table(meta$condition, meta$study))
  rio::export(sample_dist_study, file.path(PROJECT_ROOT, "results/qc/sample_distribution.csv"))
} else {
  rio::export(sample_dist, file.path(PROJECT_ROOT, "results/qc/sample_distribution.csv"))
}

log_entry("QC_SAMPLES", "OK", paste("N total:", nrow(meta)))

# 9.10 Salvar QC summary
sink(file.path(PROJECT_ROOT, "results/qc/qc_summary.md"))
cat("# QC Summary — LIHC Liver Analysis\n\n")
cat(sprintf("Total amostras fígado: %d\n", nrow(meta)))
cat(sprintf("  Normal (GTEx): %d\n", sum(meta$condition == "Normal")))
cat(sprintf("  LIHC (TCGA): %d\n", sum(meta$condition == "LIHC")))
cat(sprintf("Total genes: %d\n", ncol(expr)))
cat(sprintf("Valores NA: %d (%.2f%%)\n", na_count, na_pct))
cat(sprintf("Genes zero-var removidos: %d\n", length(zero_var_genes)))
cat(sprintf("Genes duplicados: %d\n", qc_report$dup_genes))
cat(sprintf("Range expressão: %.2f a %.2f\n", expr_range[1], expr_range[2]))
cat(sprintf("Dados integer-like: %s\n", is_integer_like))
cat(sprintf("Dados log-scale: %s\n", is_log_scale))
cat(sprintf("Batch TCGA/GTEx: %s\n",
            if ("study" %in% colnames(meta)) "2 estudos" else "Não identificado"))
cat(sprintf("Decisão: %s\n",
            if (is_integer_like || expr_range[2] > 50) "voom + limma" else "limma direto"))
sink()

log_entry("QC_DONE", "OK", "QC pré-análise concluído")

# ------------------------------------------------------------------
# 10) DECISÃO METODOLÓGICA: limma direto vs edgeR + voom + limma
# ------------------------------------------------------------------
# Se dados parecem contagens inteiras (não log), usar voom
# Se dados estão em escala log2 já normalizada, usar limma direto
use_voom <- is_integer_like && expr_range[2] > 25

log_entry("METHOD_DECISION", "OK",
          paste("use_voom:", use_voom,
                if (use_voom) "edgeR + voom + limma" else "limma direto em log2"))

# ------------------------------------------------------------------
# 11) OBC TER GENES DA VIA KEGG hsa05417
# ------------------------------------------------------------------
cat("\n--- Obtendo genes da via KEGG hsa05417 ---\n")
pathway_id <- "hsa05417"

tryCatch({
  pw <- keggGet(pathway_id)
  stopifnot(length(pw) == 1L)
  
  pathway_genes_raw <- pw[[1]]$GENE
  if (is.null(pathway_genes_raw) || length(pathway_genes_raw) == 0L) {
    stop("Nenhum gene retornado pela KEGG para ", pathway_id)
  }
  
  gene_symbols_kegg <- pathway_genes_raw[seq(2, length(pathway_genes_raw), by = 2)] |>
    str_remove("\\s*\\[.*$") |>
    str_trim() |>
    str_extract("^[^;]+") |>
    str_trim() |>
    unique() |>
    sort()
  
  genes_via_tbl <- tibble(
    pathway_id  = pathway_id,
    pathway_name = "Lipid and Atherosclerosis",
    gene_symbol = gene_symbols_kegg
  )
  
  genes_via <- genes_via_tbl$gene_symbol
  log_entry("KEGG_GENES", "OK", paste("Total via LA:", length(genes_via)))
  
  rio::export(genes_via_tbl, file.path(PROJECT_ROOT, "results/tables/KEGG_hsa05417_genes.csv"))
  
}, error = function(e) {
  log_entry("KEGG_GENES", "FAIL", paste("Erro KEGG:", e$message))
  stop("Falha ao obter genes KEGG: ", e$message)
})

# Genes da via presentes na matriz de expressão
genes_via_in_expr <- base::intersect(genes_via, colnames(expr))
log_entry("KEGG_GENES_IN_DATA", "OK", paste("Presentes na matriz:", length(genes_via_in_expr)))

if (length(genes_via_in_expr) < 2) {
  log_entry("KEGG_GENES_IN_DATA", "FAIL",
            "Menos de 2 genes da via encontrados na matriz de expressão")
  stop("Poucos genes da via LA na matriz de expressão.")
}

cat(sprintf("Genes da via LA presentes: %d / %d\n", length(genes_via_in_expr), length(genes_via)))

# ------------------------------------------------------------------
# 12) ANÁLISE DE EXPRESSÃO DIFERENCIAL
# ------------------------------------------------------------------
log_entry("DE_START", "OK", "Iniciando análise diferencial")

# Preparar matriz: genes (linhas) x amostras (colunas)
E_all <- t(expr)
storage.mode(E_all) <- "numeric"

# Remover genes com variância zero
E_all <- E_all[apply(E_all, 1, var, na.rm = TRUE) > 0, , drop = FALSE]

# Design matrix
design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)

cat(sprintf("Matriz E: %d genes x %d amostras\n", nrow(E_all), ncol(E_all)))
cat(sprintf("Design: %d amostras x %d coeficientes\n", nrow(design), ncol(design)))
stopifnot(ncol(E_all) == nrow(design))

# Aplicar voom se necessário
if (use_voom) {
  cat("Usando edgeR + voom + limma (dados tipo contagem)\n")
  dge <- DGEList(counts = E_all, group = meta$condition)
  # Filtrar genes com baixa expressão
  keep_genes <- filterByExpr(dge, group = meta$condition, min.count = 10)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
} else {
  cat("Usando limma direto (dados em escala log2 normalizada)\n")
  fit <- lmFit(E_all, design)
}

# Contraste LIHC vs Normal
contrast_matrix <- makeContrasts(LIHC_vs_Normal = LIHC - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = use_voom, robust = TRUE)

# Extrair todos os resultados
deg <- topTable(fit2, coef = "LIHC_vs_Normal", number = Inf,
                adjust.method = "BH", sort.by = "P") |>
  tibble::rownames_to_column(var = "gene_symbol")

# Classificação
deg$regulation <- dplyr::case_when(
  deg$adj.P.Val < 0.05 & deg$logFC >  1 ~ "Up_LIHC",
  deg$adj.P.Val < 0.05 & deg$logFC < -1 ~ "Down_LIHC",
  TRUE ~ "NS"
)

cat("\n--- Resumo DEGs ---\n")
print(table(deg$regulation))

n_up   <- sum(deg$regulation == "Up_LIHC")
n_down <- sum(deg$regulation == "Down_LIHC")
n_ns   <- sum(deg$regulation == "NS")

log_entry("DE_RESULTS", "OK",
          sprintf("Up: %d | Down: %d | NS: %d | Total: %d",
                  n_up, n_down, n_ns, nrow(deg)))

# Exportar tabela completa
rio::export(deg, file.path(PROJECT_ROOT, "results/deg/DEG_LIHC_vs_Normal_full.csv"))

# DEGs da via LA
deg_la <- deg |> dplyr::filter(gene_symbol %in% genes_via)
rio::export(deg_la, file.path(PROJECT_ROOT, "results/deg/DEG_LA_pathway_only.csv"))

# DEGs significativos da via LA
deg_la_sig <- deg_la |> dplyr::filter(regulation != "NS")
rio::export(deg_la_sig, file.path(PROJECT_ROOT, "results/deg/DEG_significant_LA.csv"))

# Top up/down genes
top_up <- deg |> dplyr::filter(regulation == "Up_LIHC") |>
  dplyr::arrange(adj.P.Val) |> dplyr::slice_head(n = 50)
top_down <- deg |> dplyr::filter(regulation == "Down_LIHC") |>
  dplyr::arrange(adj.P.Val) |> dplyr::slice_head(n = 50)

rio::export(top_up, file.path(PROJECT_ROOT, "results/deg/top_up_genes.csv"))
rio::export(top_down, file.path(PROJECT_ROOT, "results/deg/top_down_genes.csv"))

log_entry("DE_EXPORT", "OK", "Tabelas DEG exportadas")

# ------------------------------------------------------------------
# 13) VOLCANO PLOT — APENAS GENES DA VIA LA
# ------------------------------------------------------------------
log_entry("VOLCANO_START", "OK", "Gerando volcano plot da via LA")

lfc_cutoff <- 1
fdr_cutoff <- 0.05

deg_la$volcano_class <- dplyr::case_when(
  deg_la$adj.P.Val < fdr_cutoff & deg_la$logFC >  lfc_cutoff ~ "Up_LIHC",
  deg_la$adj.P.Val < fdr_cutoff & deg_la$logFC < -lfc_cutoff ~ "Down_LIHC",
  TRUE ~ "NS"
)

cat(sprintf("Volcano LA: Up=%d Down=%d NS=%d\n",
            sum(deg_la$volcano_class == "Up_LIHC"),
            sum(deg_la$volcano_class == "Down_LIHC"),
            sum(deg_la$volcano_class == "NS")))

# Rotular apenas genes significativos — TOP 10 up + TOP 10 down por |logFC|
genes_to_label <- deg_la |>
  dplyr::filter(volcano_class != "NS") |>
  dplyr::group_by(volcano_class) |>
  dplyr::arrange(dplyr::desc(abs(logFC))) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup() |>
  dplyr::pull(gene_symbol)

deg_la$label <- ifelse(deg_la$gene_symbol %in% genes_to_label,
                        deg_la$gene_symbol, "")

p_volcano <- ggplot(deg_la, aes(logFC, -log10(adj.P.Val), color = volcano_class)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(
    values = c(Up_LIHC = "#D73027", Down_LIHC = "#4575B4", NS = "grey75"),
    labels = c(Up_LIHC = "Upregulated LIHC", Down_LIHC = "Downregulated LIHC", NS = "NS")
  ) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_text_repel(
    data = subset(deg_la, label != ""),
    aes(label = label),
    size = 3.5, max.overlaps = 25,
    box.padding = 0.5, point.padding = 0.3,
    segment.color = "grey50", segment.size = 0.3,
    fontface = "italic"
  ) +
  labs(
    title = "Volcano Plot — Genes da Via Lipid and Atherosclerosis (hsa05417)",
    subtitle = paste("LIHC vs Normal | FDR < 0.05 | |log2FC| > 1 |",
                     length(genes_via_in_expr), "genes da via analisados"),
    x = "log2 Fold Change (LIHC / Normal)",
    y = expression(-log[10](FDR)),
    color = "Regulação"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "bottom"
  )

print(p_volcano)

ggsave(file.path(PROJECT_ROOT, "results/volcano/Volcano_LIHC_LA_pathway.png"),
       p_volcano, width = 10, height = 8, dpi = 300)
ggsave(file.path(PROJECT_ROOT, "results/volcano/Volcano_LIHC_LA_pathway.pdf"),
       p_volcano, width = 10, height = 8, device = "pdf")

# Exportar tabela dos genes rotulados
rio::export(
  deg_la |> dplyr::filter(label != ""),
  file.path(PROJECT_ROOT, "results/volcano/Volcano_labeled_genes.csv")
)

log_entry("VOLCANO_DONE", "OK", "Volcano plot exportado")

# ------------------------------------------------------------------
# 14) REDE PPI — STRINGdb (REAL)
# ------------------------------------------------------------------
log_entry("PPI_START", "OK", "Construindo rede PPI via STRINGdb")

# 14.1 Identificar DEGs significativos da via LA
deg_la_up   <- deg_la_sig |> dplyr::filter(regulation == "Up_LIHC") |> dplyr::pull(gene_symbol)
deg_la_down <- deg_la_sig |> dplyr::filter(regulation == "Down_LIHC") |> dplyr::pull(gene_symbol)
deg_la_all  <- unique(c(deg_la_up, deg_la_down))

cat(sprintf("Genes DEG LA significativos: %d (Up: %d, Down: %d)\n",
            length(deg_la_all), length(deg_la_up), length(deg_la_down)))

if (length(deg_la_all) < 3) {
  log_entry("PPI", "SKIP", "Menos de 3 DEGs LA significativos. PPI não construída.")
  
  # Criar arquivos vazios
  write.csv(data.frame(), file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_edges_score700.csv"))
  write.csv(data.frame(), file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_nodes_score700.csv"))
  
} else {
  
  # 14.2 Mapear genes para STRING IDs (com retry)
  tryCatch({
    # Aumentar timeout para downloads do STRING
    old_timeout <- getOption("timeout")
    options(timeout = 600)
    
    # Inicializar STRINGdb
    string_db <- STRINGdb$new(
      version = "12.0",
      species = 9606,           # Homo sapiens
      score_threshold = 700,
      network_type = "physical" # Interações físicas
    )
    
    log_entry("PPI_STRINGDB", "OK", "STRINGdb v12.0 inicializado | score >= 700 | physical")
    
    # Mapear símbolos para STRING IDs
    deg_string_df <- data.frame(gene_symbol = deg_la_all, stringsAsFactors = FALSE)
    
    mapped <- tryCatch({
      string_db$map(deg_string_df, "gene_symbol", removeUnmappedRows = TRUE)
    }, error = function(e) {
      log_entry("PPI_MAP", "WARN", paste("Erro no map direto:", e$message))
      NULL
    })
    
    if (is.null(mapped) || nrow(mapped) == 0) {
      mapped <- string_db$map(deg_string_df, "gene_symbol", removeUnmappedRows = FALSE)
    }
    
    mapped_ids <- mapped$STRING_id[!is.na(mapped$STRING_id)]
    
    # Genes não mapeados
    unmapped <- setdiff(deg_la_all, mapped$gene_symbol[!is.na(mapped$STRING_id)])
    if (length(unmapped) > 0) {
      unmapped_df <- data.frame(
        gene_symbol = unmapped,
        reason = "Não encontrado no STRING v12.0",
        stringsAsFactors = FALSE
      )
      rio::export(unmapped_df, file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_unmapped_genes.csv"))
      log_entry("PPI_UNMAPPED", "INFO", paste(length(unmapped), "genes não mapeados no STRING"))
    }
    
    cat(sprintf("Genes mapeados no STRING: %d / %d\n", length(mapped_ids), length(deg_la_all)))
    
    # 14.3 Obter rede de interações — download MANUAL do STRING (robusto)
    if (length(mapped_ids) >= 2) {
      
      # Baixar o arquivo completo de links físicos (score >= 700, v12.0)
      links_url <- "https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz"
      links_tmp <- tempfile(fileext = ".txt.gz")
      
      cat("Baixando links físicos STRING v12.0 (score >= 700)...\n")
      dl_ok <- tryCatch({
        download.file(links_url, links_tmp, mode = "wb")
        TRUE
      }, error = function(e) {
        log_entry("PPI_DL", "FAIL", paste("Falha download links:", e$message))
        FALSE
      })
      
      if (dl_ok) {
        # Ler o arquivo de links
        links_raw <- read.table(gzfile(links_tmp), header = TRUE, 
                                stringsAsFactors = FALSE, sep = " ")
        cat(sprintf("Links totais STRING (pré-filtro): %d\n", nrow(links_raw)))
        
        # Filtrar: protein1 OU protein2 está nos nossos mapped_ids
        links_filtered <- links_raw |>
          dplyr::filter(protein1 %in% mapped_ids & protein2 %in% mapped_ids) |>
          dplyr::filter(combined_score >= 700)
        
        cat(sprintf("Interações entre nossos genes (score >= 700): %d\n", nrow(links_filtered)))
        
        # Também incluir interações com até 50 interactores adicionais (para contexto da rede)
        links_extended <- links_raw |>
          dplyr::filter(
            (protein1 %in% mapped_ids | protein2 %in% mapped_ids) &
            combined_score >= 700
          ) |>
          dplyr::slice_max(order_by = combined_score, n = 500)
        
        # Construir interações no formato esperado
        interactions_700 <- links_filtered |>
          dplyr::transmute(
            from = protein1,
            to = protein2,
            combined_score = combined_score / 1000
          )
        
        log_entry("PPI_EDGES", "OK", paste(nrow(interactions_700), 
                  "interações STRING (download manual)"))
        
        # Limpar arquivo temporário
        unlink(links_tmp)
        
      } else {
        log_entry("PPI_EDGES", "FAIL", "Download manual STRING falhou")
        interactions_700 <- data.frame(
          from = character(), to = character(), combined_score = numeric(),
          stringsAsFactors = FALSE
        )
      }
      
      # Restaurar timeout
      options(timeout = old_timeout)
      
      # 14.4 Construir tabela de nós (mapear ENSP -> gene symbol)
      # O STRING_id no mapped tem formato "9606.ENSP0000..."
      string_to_symbol <- mapped |>
        dplyr::filter(!is.na(STRING_id)) |>
        dplyr::select(STRING_id, gene_symbol)
      
      # Identificar todos os nós envolvidos nas interações
      all_ensp <- unique(c(interactions_700$from, interactions_700$to))
      
      nodes_df <- data.frame(
        STRING_id = all_ensp,
        stringsAsFactors = FALSE
      ) |>
        dplyr::left_join(string_to_symbol, by = "STRING_id") |>
        dplyr::mutate(
          gene_symbol = ifelse(is.na(gene_symbol), STRING_id, gene_symbol)
        )
      
      # Adicionar regulação
      nodes_df$regulation <- dplyr::case_when(
        nodes_df$gene_symbol %in% deg_la_up   ~ "Up_LIHC",
        nodes_df$gene_symbol %in% deg_la_down ~ "Down_LIHC",
        TRUE ~ "Interactor"
      )
      
      # Adicionar anotação de DEG LA
      nodes_df$is_DEG_LA <- nodes_df$gene_symbol %in% deg_la_all
      
      # Converter arestas ENSP -> gene symbols
      edges_df <- interactions_700 |>
        dplyr::left_join(
          string_to_symbol |> dplyr::rename(from_gene = gene_symbol),
          by = c("from" = "STRING_id")
        ) |>
        dplyr::left_join(
          string_to_symbol |> dplyr::rename(to_gene = gene_symbol),
          by = c("to" = "STRING_id")
        ) |>
        dplyr::mutate(
          from_gene = ifelse(is.na(from_gene), from, from_gene),
          to_gene   = ifelse(is.na(to_gene), to, to_gene)
        )
      
      # 14.5 Construir grafo igraph
      g_ppi <- igraph::graph_from_data_frame(
        d = edges_df |> dplyr::transmute(from = from_gene, to = to_gene,
                                          weight = combined_score / 1000),
        directed = FALSE,
        vertices = nodes_df |> dplyr::transmute(
          name = gene_symbol, regulation, is_DEG_LA
        )
      )
      
      # Remover laços e arestas múltiplas
      g_ppi <- igraph::simplify(g_ppi, remove.loops = TRUE,
                                 edge.attr.comb = list(weight = "max"))
      
      # 14.6 Topologia
      topo <- data.frame(
        gene_symbol = igraph::V(g_ppi)$name,
        degree      = igraph::degree(g_ppi),
        betweenness = round(igraph::betweenness(g_ppi, normalized = TRUE), 4),
        closeness   = round(igraph::closeness(g_ppi, normalized = TRUE), 4),
        eigenvector = round(igraph::eigen_centrality(g_ppi)$vector, 4),
        regulation  = igraph::V(g_ppi)$regulation,
        is_DEG_LA   = igraph::V(g_ppi)$is_DEG_LA,
        stringsAsFactors = FALSE
      ) |>
        dplyr::arrange(dplyr::desc(degree))
      
      rio::export(topo, file.path(PROJECT_ROOT, "results/ppi/PPI_topology_metrics.csv"))
      
      # Identificar hubs (top 20% por degree)
      hub_threshold <- quantile(topo$degree, 0.80)
      topo$is_hub <- topo$degree >= hub_threshold
      
      cat(sprintf("Hubs identificados: %d (degree >= %.1f)\n",
                  sum(topo$is_hub), hub_threshold))
      log_entry("PPI_HUBS", "OK", paste(sum(topo$is_hub), "hubs (degree >=", 
                                        round(hub_threshold, 1), ")"))
      
      # 14.7 Componente conexo
      comps <- igraph::components(g_ppi)
      giant_idx <- which.max(comps$csize)
      cat(sprintf("Componentes: %d | Maior: %d nós | Isolados: %d\n",
                  comps$no, comps$csize[giant_idx],
                  sum(comps$csize == 1)))
      
      # Salvar tabelas PPI ANTES do plot (para garantir que sejam salvas)
      rio::export(edges_df, file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_edges_score700.csv"))
      rio::export(nodes_df, file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_nodes_score700.csv"))
      
      # Extrair maior componente para visualização
      g_giant <- igraph::induced_subgraph(g_ppi,
                                           igraph::V(g_ppi)[comps$membership == giant_idx])
      
      # 14.8 Plot da rede PPI (pré-computar atributos para evitar erros de aesthetic)
      igraph::V(g_giant)$node_degree <- igraph::degree(g_giant)
      igraph::V(g_giant)$node_label <- ifelse(
        igraph::V(g_giant)$node_degree >= quantile(igraph::V(g_giant)$node_degree, 0.70),
        igraph::V(g_giant)$name, ""
      )
      
      set.seed(4721)
      
      p_ppi <- ggraph(g_giant, layout = "stress") +
        geom_edge_link(aes(width = weight, alpha = weight), color = "grey50") +
        scale_edge_width(range = c(0.2, 2)) +
        scale_edge_alpha(range = c(0.15, 0.6)) +
        geom_node_point(aes(color = regulation, size = node_degree),
                        alpha = 0.9) +
        scale_color_manual(
          values = c(Up_LIHC = "#D73027", Down_LIHC = "#4575B4", Interactor = "grey70"),
          name = "Regulação"
        ) +
        scale_size(range = c(3, 10), name = "Degree") +
        geom_node_text(
          aes(label = node_label),
          repel = TRUE, size = 3.5, fontface = "italic",
          bg.color = "white", bg.r = 0.1
        ) +
        labs(
          title = "Rede PPI — STRING v12.0 (score ≥ 700)",
          subtitle = paste("Genes DEG da via Lipid and Atherosclerosis em LIHC |",
                           igraph::vcount(g_giant), "nós |", igraph::ecount(g_giant), "arestas"),
          caption = "Interações físicas | Fonte: STRINGdb v12.0 | Homo sapiens"
        ) +
        theme_void(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 10, color = "grey40"),
          plot.caption = element_text(size = 8, color = "grey60"),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
        )
      
      print(p_ppi)
      
      ggsave(file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_network.png"),
             p_ppi, width = 12, height = 10, dpi = 300)
      
      # 14.9 Exportar tabelas PPI
      rio::export(edges_df, file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_edges_score700.csv"))
      rio::export(nodes_df, file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_nodes_score700.csv"))
      
      # 14.10 Tabela de anotação funcional STRING dos top 10 genes (5 up + 5 down)
      cat("\n--- Anotações funcionais STRING para top 10 genes ---\n")
      top5_up <- deg_la_sig |>
        dplyr::filter(regulation == "Up_LIHC") |>
        dplyr::arrange(adj.P.Val) |>
        dplyr::slice_head(n = 5)
      top5_down <- deg_la_sig |>
        dplyr::filter(regulation == "Down_LIHC") |>
        dplyr::arrange(adj.P.Val) |>
        dplyr::slice_head(n = 5)
      top10_genes <- dplyr::bind_rows(top5_up, top5_down)
      
      # Fallback manual com anotações funcionais curadas do STRING/UniProt
      known <- data.frame(
        gene_symbol = c("BAX","HSP90AB1","MMP1","MMP9","PPARG","CD36","ABCG1",
                        "IKBKE","LY96","VCAM1","CALML5","CXCL2","CAMK2B","CALML6",
                        "MIB2","IKBKB","STAT3","ABCA1","TNFRSF10B","SELE","CALML3",
                        "IRF7","FOS","TLR2","CYP2C8","CYP2B6","CXCL3","CYP2A7",
                        "APOA4","CYP2C9","HSPA6","HSPA1B","MAPK10","MAP3K5","SELP",
                        "PIK3R2","LBP","CYP2A6","CYP1A1","PYCARD","PLCB2","CAMK2A"),
        func = c(
          "Regulador de apoptose; permeabilização da membrana mitocondrial; família BCL-2",
          "Chaperona molecular; dobramento proteico; família HSP90; estabilidade proteica",
          "Colagenase; degradação da matriz extracelular; metaloproteinase; remodelamento tecidual",
          "Gelatinase; degradação da matriz extracelular; metaloproteinase; invasão tumoral",
          "Receptor nuclear; metabolismo lipídico; diferenciação de adipócitos; anti-inflamatório",
          "Receptor scavenger; transporte de ácidos graxos; captação lipídica; adesão celular",
          "Efluxo de colesterol; transporte lipídico; transportador ABC; biogênese de HDL",
          "Quinase; sinalização NF-kappaB; resposta antiviral; via TBK1/IKKi",
          "Correceptor TLR4; reconhecimento de LPS; imunidade inata; proteína MD-2",
          "Adesão celular; migração leucocitária; ativação endotelial; superfamília Ig",
          "Sinalização de cálcio; calmodulina-like; domínio EF-hand; proliferação celular",
          "Quimiocina; resposta inflamatória; quimiotaxia de neutrófilos; resposta imune",
          "Quinase dependente de Ca²⁺/calmodulina; plasticidade sináptica; família CaMK",
          "Sinalização de cálcio; calmodulina-like; domínio EF-hand",
          "E3 ubiquitina ligase; sinalização Notch; ubiquitinação proteica",
          "Sinalização NF-kappaB; complexo IKK; resposta inflamatória; regulação imune",
          "Fator de transcrição; via JAK-STAT; proliferação celular; regulação imune",
          "Efluxo de colesterol; biogênese de HDL; transporte reverso de colesterol",
          "Receptor de morte; superfamília TNF; apoptose extrínseca; receptor TRAIL",
          "Adesão celular; migração leucocitária; ativação endotelial; selectina",
          "Sinalização de cálcio; calmodulina-like; domínio EF-hand",
          "Fator de transcrição; sinalização de interferon; resposta antiviral; imunidade inata",
          "Fator de transcrição; complexo AP-1; proto-oncogene; resposta ao estresse",
          "Receptor Toll-like; imunidade inata; reconhecimento de padrões; resposta inflamatória",
          "Citocromo P450; metabolismo xenobiótico; metabolismo do ácido araquidônico",
          "Citocromo P450; metabolismo xenobiótico; metabolismo de drogas/esteroides",
          "Quimiocina; resposta inflamatória; ativação de neutrófilos; ligante CXCR2",
          "Citocromo P450; metabolismo xenobiótico; metabolismo da cumarina",
          "Apolipoproteína; transporte lipídico; componente HDL; metabolismo do colesterol",
          "Citocromo P450; metabolismo xenobiótico; metabolismo de drogas",
          "Proteína de choque térmico; chaperona; resposta ao estresse; família HSP70",
          "Proteína de choque térmico; chaperona; resposta ao estresse; família HSP70",
          "MAP quinase; sinalização JNK; resposta ao estresse; apoptose neuronal",
          "MAP quinase quinase quinase; ASK1; sinalização de estresse; via JNK/p38",
          "Adesão celular; migração leucocitária; ativação plaquetária; selectina",
          "Subunidade regulatória PI3K; sinalização de insulina; crescimento celular",
          "Proteína ligante de LPS; imunidade inata; ativação de TLR4; fase aguda",
          "Citocromo P450; metabolismo xenobiótico; metabolismo da cumarina",
          "Citocromo P450; metabolismo xenobiótico; metabolismo de hidrocarbonetos arílicos",
          "Componente do inflamassomo; apoptose; ativação de caspases; imunidade inata",
          "Fosfolipase C; sinalização de cálcio; transdução de sinal intracelular",
          "Quinase dependente de Ca²⁺/calmodulina; plasticidade sináptica; família CaMK"
        ), stringsAsFactors = FALSE)
      top10_func <- top10_genes |>
        dplyr::left_join(known, by = "gene_symbol") |>
        dplyr::rename(STRING_function = func) |>
        dplyr::select(gene_symbol, logFC, adj.P.Val, regulation, STRING_function)
      
      rio::export(top10_func, file.path(PROJECT_ROOT, "results/ppi/Top10_genes_STRING_functions.csv"))
      log_entry("PPI_ANNOT", "OK", paste(nrow(top10_func), "genes anotados com funções STRING"))
      
      # 14.11 Registro de consistência da PPI
      ppi_log_df <- data.frame(
        timestamp = as.character(Sys.time()),
        n_input = length(deg_la_all),
        n_mapped = length(mapped_ids),
        n_interact = nrow(interactions_700),
        n_nodes = igraph::vcount(g_ppi),
        n_edges = igraph::ecount(g_ppi),
        n_components = comps$no,
        largest_comp = comps$csize[giant_idx],
        hubs_n = sum(topo$is_hub, na.rm = TRUE),
        hubs = paste(topo$gene_symbol[topo$is_hub], collapse = "; "),
        stringsAsFactors = FALSE
      )
      rio::export(ppi_log_df, file.path(PROJECT_ROOT, "results/ppi/PPI_consistency_log.csv"))
      log_entry("PPI_CONSISTENCY", "OK", paste("Hubs:", ppi_log_df$hubs))
      
      log_entry("PPI_DONE", "OK",
                paste("Rede:", igraph::vcount(g_ppi), "nós |", igraph::ecount(g_ppi),
                      "arestas | Maior comp:", igraph::vcount(g_giant), "nós"))
      
    } else {
      log_entry("PPI", "SKIP", "Menos de 2 genes mapeados no STRING. Rede não construída.")
      write.csv(data.frame(), file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_edges_score700.csv"))
      write.csv(data.frame(), file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_nodes_score700.csv"))
    }
    
  }, error = function(e) {
    log_entry("PPI_ERROR", "FAIL", paste("Erro STRINGdb:", e$message))
    cat("AVISO: STRINGdb falhou. Erro:", e$message, "\n")
    write.csv(data.frame(error = e$message),
              file.path(PROJECT_ROOT, "results/ppi/PPI_STRING_error.csv"))
  })
}

# ------------------------------------------------------------------
# 15) GSEA (rank-based) + GSVA (por amostra) — PRIMEIROS PROCEDIMENTOS
# ------------------------------------------------------------------
map_to_entrez <- function(symbols) {
  tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE) |>
      dplyr::pull(ENTREZID)
  }, error = function(e) {
    log_entry("ENRICH_MAP", "WARN", paste("Erro mapeamento ENTREZ:", e$message))
    return(character(0))
  })
}

log_entry("GSEA_START", "OK", "GSEA rank-based + GSVA por amostra")

# Lista ranqueada por logFC (todos os genes testados)
ranked_deg <- deg |>
  dplyr::filter(!is.na(logFC)) |>
  dplyr::arrange(dplyr::desc(logFC))

# Vetor ranqueado por símbolo (para fgsea) e por ENTREZ (para gseKEGG)
sym_rank_vec <- ranked_deg$logFC
names(sym_rank_vec) <- ranked_deg$gene_symbol
sym_rank_vec <- sym_rank_vec[!duplicated(names(sym_rank_vec))]

entrez_map_df <- tryCatch(
  bitr(ranked_deg$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db, drop = TRUE),
  error = function(e) data.frame(SYMBOL = character(), ENTREZID = character()))
if (nrow(entrez_map_df) > 0) {
  entrez_df <- data.frame(SYMBOL = ranked_deg$gene_symbol,
                          logFC  = ranked_deg$logFC,
                          stringsAsFactors = FALSE) |>
    dplyr::inner_join(entrez_map_df, by = "SYMBOL") |>
    dplyr::distinct(ENTREZID, .keep_all = TRUE) |>
    dplyr::arrange(dplyr::desc(logFC))
  entrez_rank_vec <- entrez_df$logFC
  names(entrez_rank_vec) <- entrez_df$ENTREZID
} else {
  entrez_rank_vec <- NULL
}

# GSEA KEGG
if (!is.null(entrez_rank_vec) && length(entrez_rank_vec) >= 10) {
  gsea_kegg <- tryCatch(
    gseKEGG(geneList = entrez_rank_vec, organism = "hsa",
            eps = 0, pvalueCutoff = 0.1, seed = 4721),
    error = function(e) NULL)
  if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
    gsea_kegg_df <- as.data.frame(gsea_kegg)
    gsea_kegg_df <- gsea_kegg_df[!is.na(gsea_kegg_df$ID), , drop = FALSE]
    rio::export(gsea_kegg_df, file.path(PROJECT_ROOT, "results/enrichment/GSEA_KEGG.csv"))
    log_entry("GSEA_KEGG", "OK", paste(nrow(gsea_kegg_df), "conjuntos testados"))
  } else {
    rio::export(data.frame(), file.path(PROJECT_ROOT, "results/enrichment/GSEA_KEGG.csv"))
    log_entry("GSEA_KEGG", "INFO", "Nenhum conjunto enriquecido")
  }
}

# GSEA Hallmark (MSigDB) via fgsea
if (length(sym_rank_vec) >= 10) {
  hall <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |>
    dplyr::select(gs_name, gene_symbol)
  hall_list <- split(hall$gene_symbol, hall$gs_name)
  
  set.seed(4721)
  fgsea_res <- fgsea::fgsea(pathways = hall_list, stats = sym_rank_vec,
                            minSize = 5, maxSize = 500, nPermSimple = 10000)
  if (nrow(fgsea_res) > 0) {
    fgsea_out <- fgsea_res |>
      dplyr::mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";",
                                         FUN.VALUE = character(1))) |>
      dplyr::arrange(padj)
    rio::export(fgsea_out, file.path(PROJECT_ROOT, "results/enrichment/GSEA_HALLMARK.csv"))
    log_entry("GSEA_HALLMARK", "OK", paste(nrow(fgsea_out), "conjuntos testados"))
  } else {
    rio::export(data.frame(), file.path(PROJECT_ROOT, "results/enrichment/GSEA_HALLMARK.csv"))
    log_entry("GSEA_HALLMARK", "INFO", "Nenhum conjunto testado")
  }
}

# GSVA — escore de enriquecimento por amostra
if (exists("hall_list") && length(hall_list) > 0) {
  gsva_sets <- hall_list[grep("^HALLMARK", names(hall_list))]
} else {
  gsva_sets <- list()
}
gsva_sets[["KEGG_hsa05417_Lipid_Atherosclerosis"]] <- genes_via_in_expr
gsva_sets <- lapply(gsva_sets, function(g) base::intersect(g, rownames(E_all)))
gsva_sets <- gsva_sets[lengths(gsva_sets) >= 3]

if (length(gsva_sets) >= 2) {
  gsva_res <- tryCatch({
    param <- GSVA::gsvaParam(exprData = E_all, geneSets = gsva_sets,
                             minSize = 3, maxSize = 500)
    GSVA::gsva(param)
  }, error = function(e) NULL)
  
  if (!is.null(gsva_res) && nrow(gsva_res) > 0) {
    gsva_mat <- as.data.frame(gsva_res)
    rio::export(gsva_mat, file.path(PROJECT_ROOT, "results/enrichment/GSVA_scores.csv"))
    
    # Comparação LIHC vs Normal por conjunto (teste t)
    gsva_summary <- lapply(rownames(gsva_res), function(gs) {
      scores <- as.numeric(gsva_res[gs, ])
      li <- scores[meta$condition == "LIHC"]
      no <- scores[meta$condition == "Normal"]
      tt <- tryCatch(t.test(li, no), error = function(e) NULL)
      if (is.null(tt)) {
        c(pathway = gs, diff_score = NA_real_, p.value = NA_real_)
      } else {
        c(pathway = gs, diff_score = unname(tt$estimate[1] - tt$estimate[2]),
          p.value = tt$p.value)
      }
    })
    gsva_summary <- as.data.frame(do.call(rbind, gsva_summary), stringsAsFactors = FALSE)
    gsva_summary$diff_score <- as.numeric(gsva_summary$diff_score)
    gsva_summary$p.value <- as.numeric(gsva_summary$p.value)
    gsva_summary$padj <- p.adjust(gsva_summary$p.value, method = "BH")
    gsva_summary <- gsva_summary |> dplyr::arrange(p.value)
    rio::export(gsva_summary, file.path(PROJECT_ROOT, "results/enrichment/GSVA_summary.csv"))
    log_entry("GSVA", "OK", paste(nrow(gsva_summary), "conjuntos avaliados"))
  } else {
    rio::export(data.frame(), file.path(PROJECT_ROOT, "results/enrichment/GSVA_scores.csv"))
    log_entry("GSVA", "INFO", "GSVA não produziu resultados")
  }
} else {
  log_entry("GSVA", "SKIP", "Conjuntos gênicos insuficientes para GSVA")
}

log_entry("GSEA_DONE", "OK", "GSEA + GSVA concluídos")

# ------------------------------------------------------------------
# 16) ENRIQUECIMENTO FUNCIONAL (ORA — analise complementar)
# ------------------------------------------------------------------
log_entry("ENRICH_START", "OK", "Enriquecimento KEGG/GO (análise complementar)")

# Background: genoma completo (ORA padrão — corrige o resultado vazio da
# versão anterior, que restringia o universo aos 212 genes da via)

# Up regulated
entrez_up <- map_to_entrez(top_up$gene_symbol)
entrez_down <- map_to_entrez(top_down$gene_symbol)

cat(sprintf("ENTREZ mapeados: Up=%d Down=%d\n",
            length(entrez_up), length(entrez_down)))

# KEGG enrichment
if (length(entrez_up) >= 5) {
  kegg_up <- enrichKEGG(gene = entrez_up, organism = "hsa",
                         pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
    rio::export(as.data.frame(kegg_up), file.path(PROJECT_ROOT, "results/enrichment/KEGG_up.csv"))
    log_entry("ENRICH_KEGG_UP", "OK", paste(nrow(kegg_up@result), "vias"))
  } else {
    log_entry("ENRICH_KEGG_UP", "INFO", "Nenhuma via significativa")
    rio::export(data.frame(), file.path(PROJECT_ROOT, "results/enrichment/KEGG_up.csv"))
  }
} else {
  log_entry("ENRICH_KEGG_UP", "SKIP", "Poucos genes up para enriquecimento")
}

if (length(entrez_down) >= 5) {
  kegg_down <- enrichKEGG(gene = entrez_down, organism = "hsa",
                           pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
    rio::export(as.data.frame(kegg_down), file.path(PROJECT_ROOT, "results/enrichment/KEGG_down.csv"))
    log_entry("ENRICH_KEGG_DOWN", "OK", paste(nrow(kegg_down@result), "vias"))
  } else {
    log_entry("ENRICH_KEGG_DOWN", "INFO", "Nenhuma via significativa")
    rio::export(data.frame(), file.path(PROJECT_ROOT, "results/enrichment/KEGG_down.csv"))
  }
} else {
  log_entry("ENRICH_KEGG_DOWN", "SKIP", "Poucos genes down para enriquecimento")
}

# GO Biological Process enrichment
if (length(entrez_up) >= 5) {
  go_up <- enrichGO(gene = entrez_up, OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE)
  if (!is.null(go_up) && nrow(go_up@result) > 0) {
    rio::export(as.data.frame(go_up), file.path(PROJECT_ROOT, "results/enrichment/GO_BP_up.csv"))
    log_entry("ENRICH_GO_UP", "OK", paste(nrow(go_up@result), "termos GO"))
  } else {
    rio::export(data.frame(), file.path(PROJECT_ROOT, "results/enrichment/GO_BP_up.csv"))
    log_entry("ENRICH_GO_UP", "INFO", "Nenhum termo GO significativo")
  }
}

if (length(entrez_down) >= 5) {
  go_down <- enrichGO(gene = entrez_down, OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                      readable = TRUE)
  if (!is.null(go_down) && nrow(go_down@result) > 0) {
    rio::export(as.data.frame(go_down), file.path(PROJECT_ROOT, "results/enrichment/GO_BP_down.csv"))
    log_entry("ENRICH_GO_DOWN", "OK", paste(nrow(go_down@result), "termos GO"))
  } else {
    rio::export(data.frame(), file.path(PROJECT_ROOT, "results/enrichment/GO_BP_down.csv"))
    log_entry("ENRICH_GO_DOWN", "INFO", "Nenhum termo GO significativo")
  }
}

# Dotplot combinado (se houver resultados)
kegg_up_df <- tryCatch(as.data.frame(kegg_up), error = function(e) NULL)
kegg_down_df <- tryCatch(as.data.frame(kegg_down), error = function(e) NULL)

if (!is.null(kegg_up_df) && nrow(kegg_up_df) > 0) kegg_up_df$set <- "Up_LIHC"
if (!is.null(kegg_down_df) && nrow(kegg_down_df) > 0) kegg_down_df$set <- "Down_LIHC"

kegg_combined <- dplyr::bind_rows(kegg_up_df, kegg_down_df)

if (nrow(kegg_combined) > 0) {
  kegg_combined <- kegg_combined |>
    dplyr::filter(p.adjust < 0.05) |>
    dplyr::mutate(logFDR = -log10(p.adjust))
  
  if (nrow(kegg_combined) > 0) {
    p_enrich <- ggplot(kegg_combined,
                       aes(x = logFDR, y = reorder(Description, logFDR),
                           size = Count, color = logFDR)) +
      geom_point(alpha = 0.9) +
      facet_wrap(~ set, scales = "free_y") +
      scale_color_gradient(low = "#74add1", high = "#d73027", name = "-log10(FDR)") +
      scale_size(range = c(3, 9), name = "Gene Count") +
      labs(x = "-log10(FDR)", y = "KEGG Pathway",
           title = "KEGG Enrichment — Genes DEG LIHC") +
      theme_bw(base_size = 13) +
      theme(axis.text.y = element_text(size = 8),
            strip.text = element_text(size = 12, face = "bold"))
    
    ggsave(file.path(PROJECT_ROOT, "results/enrichment/enrichment_dotplot.png"),
           p_enrich, width = 12, height = 8, dpi = 300)
  }
}

log_entry("ENRICH_DONE", "OK", "Enriquecimento concluído")

# ------------------------------------------------------------------
# 17) HEATMAP — TOP DEGs
# ------------------------------------------------------------------
cat("\n--- Heatmap Top DEGs ---\n")
top_deg_genes <- deg |>
  dplyr::filter(regulation != "NS") |>
  dplyr::arrange(adj.P.Val) |>
  dplyr::slice_head(n = 50) |>
  dplyr::pull(gene_symbol)

top_deg_in_expr <- base::intersect(top_deg_genes, rownames(E_all))

if (length(top_deg_in_expr) >= 5) {
  mat_heat <- E_all[top_deg_in_expr, , drop = FALSE]
  mat_scaled <- t(scale(t(mat_heat)))
  mat_scaled <- pmin(pmax(mat_scaled, -3), 3)
  
  annotation_col <- data.frame(
    Condition = meta$condition,
    row.names = colnames(mat_scaled)
  )
  
  ann_colors <- list(Condition = c(Normal = "#2E86AB", LIHC = "#A23B72"))
  
  pheatmap(mat_scaled,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           cluster_cols = TRUE,
           clustering_method = "ward.D2",
           show_rownames = TRUE,
           show_colnames = FALSE,
           fontsize_row = 6,
           main = paste("Top", length(top_deg_in_expr), "DEGs"),
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
           breaks = seq(-3, 3, length.out = 51),
           height = 10, width = 12,
           filename = file.path(PROJECT_ROOT, "results/figures/Heatmap_Top_DEGs.png"))
  
  log_entry("HEATMAP", "OK", paste(length(top_deg_in_expr), "genes"))
}

# ------------------------------------------------------------------
# 18) BENCHMARK
# ------------------------------------------------------------------
T1 <- Sys.time()
total_time <- round(as.numeric(difftime(T1, T0, units = "secs")), 2)

bench <- data.frame(
  metric = c(
    "tempo_total_execucao_seg",
    "num_amostras_total",
    "num_amostras_LIHC",
    "num_amostras_Normal",
    "num_genes_total_matriz",
    "num_genes_via_KEGG",
    "num_genes_via_na_matriz",
    "num_DEGs_up",
    "num_DEGs_down",
    "num_DEGs_total",
    "num_DEGs_LA_significativos",
    "num_genes_mapeados_STRING",
    "num_interacoes_STRING_700",
    "num_hubs_PPI",
    "metodo_DE",
    "batch_TGCA_GTEx"
  ),
  value = c(
    total_time,
    nrow(meta),
    sum(meta$condition == "LIHC"),
    sum(meta$condition == "Normal"),
    ncol(expr),
    length(genes_via),
    length(genes_via_in_expr),
    n_up,
    n_down,
    n_up + n_down,
    nrow(deg_la_sig),
    if (exists("mapped_ids")) length(mapped_ids) else 0,
    if (exists("interactions_700")) nrow(interactions_700) else 0,
    if (exists("topo")) sum(topo$is_hub, na.rm = TRUE) else 0,
    if (use_voom) "voom+limma" else "limma_direto",
    if ("study" %in% colnames(meta)) "identificado" else "nao_identificado"
  ),
  stringsAsFactors = FALSE
)

rio::export(bench, file.path(PROJECT_ROOT, "results/audit/benchmark_pipeline.csv"))

sink(file.path(PROJECT_ROOT, "results/audit/benchmark_summary.md"))
cat("# Benchmark Pipeline LIHC\n\n")
cat(sprintf("**Tempo total:** %.2f segundos\n\n", total_time))
cat(sprintf("- Amostras: %d (LIHC: %d, Normal: %d)\n",
            nrow(meta), sum(meta$condition == "LIHC"), sum(meta$condition == "Normal")))
cat(sprintf("- Genes na matriz: %d\n", ncol(expr)))
cat(sprintf("- Genes via LA (KEGG): %d\n", length(genes_via)))
cat(sprintf("- Genes via LA na matriz: %d\n", length(genes_via_in_expr)))
cat(sprintf("- DEGs (FDR<0.05, |logFC|>1): Up=%d, Down=%d, Total=%d\n", n_up, n_down, n_up + n_down))
cat(sprintf("- Método DE: %s\n", if (use_voom) "edgeR + voom + limma" else "limma direto"))
sink()

log_entry("BENCHMARK", "OK", paste("Tempo total:", total_time, "segundos"))

# ------------------------------------------------------------------
# 19) SALVAR SESSIONINFO
# ------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(PROJECT_ROOT, "results/audit/sessionInfo.txt"))

# ------------------------------------------------------------------
# 20) SALVAR WORKSPACE
# ------------------------------------------------------------------
save.image(file = file.path(PROJECT_ROOT, "results/audit/workspace_LIHC.RData"))
gc()

log_entry("FINAL", "OK", paste("Pipeline concluído com sucesso. Tempo:", total_time, "s"))

cat("\n")
cat("============================================================\n")
cat(" PIPELINE CONCLUÍDO COM SUCESSO\n")
cat("============================================================\n")
cat(sprintf("Tempo total: %.2f segundos\n", total_time))
cat(sprintf("Amostras: %d (LIHC: %d, Normal: %d)\n",
            nrow(meta), sum(meta$condition == "LIHC"), sum(meta$condition == "Normal")))
cat(sprintf("DEGs significativos: %d (Up: %d, Down: %d)\n", n_up + n_down, n_up, n_down))
cat(sprintf("Resultados em: %s\n", file.path(PROJECT_ROOT, "results")))
cat("============================================================\n")
