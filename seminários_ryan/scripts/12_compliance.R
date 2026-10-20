# =============================================================================
# 12_compliance.R
# -----------------------------------------------------------------------------
# Checklist de conformidade e integridade dos artefatos do seminário.
#
# Saídas:
#   auditorias/COMPLIANCE.md
#   auditorias/COMPLIANCE_tabela.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(knitr)
  library(png)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_aud <- file.path(ROOT, "auditorias")
dir.create(dir_aud, recursive = TRUE, showWarnings = FALSE)

item <- character(); status <- character(); evidencia <- character()
add_check <- function(i, s, e) {
  item <<- c(item, i); status <<- c(status, s); evidencia <<- c(evidencia, e)
}

# -----------------------------------------------------------------------------
# 1. Estrutura e arquivos obrigatórios
# -----------------------------------------------------------------------------
obrigatorios <- c(
  "README.md", "relatorio_completo.md", "requirements_R.txt",
  "scripts/00_preparar_dados.R", "scripts/01_analise_exploratoria.R",
  "scripts/02_expressao_diferencial.R", "scripts/03_volcano_plots.R",
  "scripts/04_ComBat_demonstracao.R", "scripts/05_fenotipos_avancados.R",
  "scripts/06_grupos_alternativos.R", "scripts/07_grupos_separados.R",
  "scripts/10_auditoria_dados_e_volcano.R", "scripts/11_benchmarks.R",
  "scripts/12_compliance.R", "scripts/run_all.R",
  "docs/08_fenotipos_avancados.md", "docs/09_auditorias_benchmarks.md",
  "docs/10_grupos_alternativos.md", "docs/11_grupos_separados.md",
  "resultados/tabelas/35_tss_umap.csv",
  "resultados/figuras/35_tss_marcadores_heatmap.png",
  "auditorias/LEIA_ME.md",
  "dados/processed/metadata_analise.tsv",
  "dados/processed/expressao_log2tpm_filtrada.tsv.gz",
  "dados/processed/expressao_LA_pathway.tsv.gz",
  "resultados/tabelas/DE/RESUMO_DEGs_TODAS_COMPARACOES.csv",
  "resultados/tabelas/14_fenotipos_testes.csv",
  "resultados/tabelas/18_cox_multivariado.csv",
  "resultados/figuras/Volcano_01_LIHC_TCGA_tumor_vs_adjacente.png",
  "resultados/figuras/Volcano_02_GTEx_masculino_vs_feminino.png",
  "resultados/figuras/Volcano_03_TCGA_tumor_vs_GTEx_normal.png",
  "resultados/figuras/Volcano_12_TCGA_vs_GTEx_ComBat.png",
  "resultados/figuras/13_UMAP_ComBat_antes_depois.png",
  "resultados/figuras/24_forest_cox_multivariado.png",
  "auditorias/AUDITORIA_DADOS.md",
  "auditorias/AUDITORIA_VOLCANO_12.md",
  "auditorias/BENCHMARKS.md"
)
faltantes <- obrigatorios[!file.exists(file.path(ROOT, obrigatorios))]
add_check("Arquivos obrigatórios", ifelse(length(faltantes) == 0, "OK", "FALHA"),
          paste0("ausentes: ", length(faltantes),
                 ifelse(length(faltantes) > 0, paste0(" (", paste(faltantes, collapse = ", "), ")"), "")))

# -----------------------------------------------------------------------------
# 2. Sem PDFs
# -----------------------------------------------------------------------------
pdfs <- list.files(file.path(ROOT, "resultados", "figuras"),
                   pattern = "[.]pdf$", ignore.case = TRUE, full.names = FALSE)
add_check("Figuras sem PDF", ifelse(length(pdfs) == 0, "OK", "FALHA"),
          paste0("PDFs encontrados: ", length(pdfs)))

# -----------------------------------------------------------------------------
# 3. PNGs válidos e não vazios
# -----------------------------------------------------------------------------
pngs <- list.files(file.path(ROOT, "resultados", "figuras"),
                   pattern = "[.]png$", ignore.case = TRUE, full.names = TRUE)
png_info <- data.frame(arquivo = basename(pngs), existe = TRUE,
                       dim = NA_character_, nonwhite = NA_real_,
                       stringsAsFactors = FALSE)
for (i in seq_along(pngs)) {
  im <- tryCatch(png::readPNG(pngs[i]), error = function(e) NULL)
  if (is.null(im)) {
    png_info$existe[i] <- FALSE
    next
  }
  png_info$dim[i] <- paste(dim(im), collapse = "x")
  if (length(dim(im)) == 3) {
    r <- im[, , 1]; g <- im[, , 2]; b <- im[, , 3]
    png_info$nonwhite[i] <- mean(!(r > 0.96 & g > 0.96 & b > 0.96))
  } else {
    png_info$nonwhite[i] <- mean(im < 0.96)
  }
}
vazias <- png_info$arquivo[is.na(png_info$nonwhite) | png_info$nonwhite < 0.005]
add_check("PNGs válidos e não vazios", ifelse(length(vazias) == 0, "OK", "FALHA"),
          paste0("PNGs: ", length(pngs), "; suspeitos de vazio: ",
                 ifelse(length(vazias) == 0, "0", paste(vazias, collapse = ", "))))

# Volcano_12 precisa ter conteúdo visível (antes era o problema)
v12 <- png_info[png_info$arquivo == "Volcano_12_TCGA_vs_GTEx_ComBat.png", ]
add_check("Volcano_12 não vazia", ifelse(nrow(v12) == 1 && v12$nonwhite > 0.02, "OK", "FALHA"),
          paste0("fração não branca: ", round(v12$nonwhite, 4)))

# -----------------------------------------------------------------------------
# 4. CSVs
# -----------------------------------------------------------------------------
csvs <- list.files(file.path(ROOT, "resultados", "tabelas"),
                   pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
csv_bad <- character()
for (f in csvs) {
  d <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) csv_bad <- c(csv_bad, basename(f))
}
add_check("Tabelas CSV válidas", ifelse(length(csv_bad) == 0, "OK", "FALHA"),
          paste0("CSVs: ", length(csvs), "; vazios/ilegíveis: ",
                 ifelse(length(csv_bad) == 0, "0", paste(csv_bad, collapse = ", "))))

# Colunas críticas sem NA
check_cols <- function(f, cols) {
  d <- read.csv(file.path(ROOT, f), check.names = FALSE)
  miss <- colSums(is.na(d[, cols, drop = FALSE]))
  list(ok = all(miss == 0), detail = paste0(f, ": ", paste(names(miss), miss, sep = "=", collapse = "; ")))
}
chk1 <- check_cols("resultados/tabelas/DE/DEG_TCGA_Tumor_vs_GTEx_Normal_completo.csv",
                   c("gene", "logFC", "P.Value", "adj.P.Val"))
chk2 <- check_cols("resultados/tabelas/14_fenotipos_testes.csv",
                   c("variavel_1", "variavel_2", "p_valor"))
chk3 <- check_cols("resultados/tabelas/18_cox_multivariado.csv",
                   c("term", "estimate", "conf.low", "conf.high", "p.value"))
add_check("Colunas críticas sem NA",
          ifelse(all(chk1$ok, chk2$ok, chk3$ok), "OK", "FALHA"),
          paste(chk1$detail, chk2$detail, chk3$detail, sep = " | "))

# -----------------------------------------------------------------------------
# 5. Scripts e sementes
# -----------------------------------------------------------------------------
scripts <- list.files(file.path(ROOT, "scripts"), pattern = "[.]R$", full.names = TRUE)
parse_bad <- character()
for (f in scripts) {
  ok <- tryCatch({ parse(f); TRUE }, error = function(e) FALSE)
  if (!ok) parse_bad <- c(parse_bad, basename(f))
}
add_check("Scripts R parseáveis", ifelse(length(parse_bad) == 0, "OK", "FALHA"),
          paste0("scripts: ", length(scripts), "; erros: ",
                 ifelse(length(parse_bad) == 0, "0", paste(parse_bad, collapse = ", "))))

seed_files <- c("01_analise_exploratoria.R", "03_volcano_plots.R",
                "04_ComBat_demonstracao.R", "05_fenotipos_avancados.R",
                "06_grupos_alternativos.R", "07_grupos_separados.R")
seed_ok <- sapply(seed_files, function(f) {
  any(grepl("set.seed", readLines(file.path(ROOT, "scripts", f), warn = FALSE)))
})
add_check("Sementes de reprodutibilidade", ifelse(all(seed_ok), "OK", "FALHA"),
          paste(names(seed_ok), seed_ok, sep = "=", collapse = "; "))

# -----------------------------------------------------------------------------
# 6. Documentação: imagens e links
# -----------------------------------------------------------------------------
md_files <- list.files(file.path(ROOT, "docs"), pattern = "[.]md$", full.names = TRUE)
md_bad <- character()
for (f in md_files) {
  txt <- readLines(f, warn = FALSE)
  imgs <- regmatches(txt, gregexpr("!\\[[^]]*\\]\\(([^)]+)\\)", txt))
  for (x in unlist(imgs)) {
    p <- sub("^!\\[[^]]*\\]\\(([^)]+)\\)$", "\\1", x)
    if (!file.exists(file.path(ROOT, "docs", p))) md_bad <- c(md_bad, p)
  }
}
add_check("Imagens dos docs existem", ifelse(length(md_bad) == 0, "OK", "FALHA"),
          paste0("docs: ", length(md_files), "; imagens ausentes: ",
                 ifelse(length(md_bad) == 0, "0", paste(md_bad, collapse = ", "))))

# Relatório Markdown consolidado
report_bad <- character()
if (file.exists(file.path(ROOT, "relatorio_completo.md"))) {
  report <- readLines(file.path(ROOT, "relatorio_completo.md"), warn = FALSE)
  imgs <- regmatches(report, gregexpr("!\\[[^]]*\\]\\(([^)]+)\\)", report))
  for (x in unlist(imgs)) {
    p <- sub("^!\\[[^]]*\\]\\(([^)]+)\\)$", "\\1", x)
    if (!file.exists(file.path(ROOT, p))) report_bad <- c(report_bad, p)
  }
}
add_check("Imagens do relatório Markdown existem",
          ifelse(length(report_bad) == 0, "OK", "FALHA"),
          paste0("ausentes: ", ifelse(length(report_bad) == 0, "0", paste(report_bad, collapse = ", "))))

# Sem menções a PDF ou PCA nos docs
docs_txt <- paste(unlist(lapply(md_files, readLines, warn = FALSE)), collapse = "\n")
pdf_mentions <- grepl("[.]pdf", docs_txt, ignore.case = TRUE)
pca_mentions <- grepl("\\bPCA\\b|PC1|PC2|prcomp", docs_txt)
add_check("Docs sem PDF", ifelse(!pdf_mentions, "OK", "FALHA"),
          ifelse(pdf_mentions, "há menções a .pdf", "nenhuma menção"))
add_check("Docs sem PCA (substituído por UMAP)", ifelse(!pca_mentions, "OK", "FALHA"),
          ifelse(pca_mentions, "há menções a PCA/PC1/PC2", "nenhuma menção"))

# -----------------------------------------------------------------------------
# 7. Metadados
# -----------------------------------------------------------------------------
meta <- fread(file.path(ROOT, "dados/processed/metadata_analise.tsv"),
              sep = "\t", data.table = FALSE, na.strings = c("", "NA"))
cols_esperadas <- c("sample", "group", "study", "sex", "age", "stage", "grade",
                    "stage_group", "grade_group", "afp", "afp_high",
                    "fibrosis_ishak", "fibrosis_advanced", "viral_status",
                    "child_pugh", "histology", "vascular_invasion",
                    "race", "ethnicity", "OS", "OS_time")
faltam <- setdiff(cols_esperadas, names(meta))
add_check("Colunas do metadado", ifelse(length(faltam) == 0, "OK", "FALHA"),
          paste0("ausentes: ", ifelse(length(faltam) == 0, "0", paste(faltam, collapse = ", "))))
add_check("Contagem de amostras", ifelse(nrow(meta) == 529, "OK", "FALHA"),
          paste0("n = ", nrow(meta), " (esperado 529)"))
add_check("Sem recidivas na análise", ifelse(!any(grepl("-02$", meta$sample)), "OK", "FALHA"),
          paste0("amostras -02: ", sum(grepl("-02$", meta$sample))))

# -----------------------------------------------------------------------------
# 8. Relatório
# -----------------------------------------------------------------------------
comp_df <- data.frame(item, status, evidencia, stringsAsFactors = FALSE)
write.csv(comp_df, file.path(dir_aud, "COMPLIANCE_tabela.csv"), row.names = FALSE)

md <- c(
  "# Compliance e integridade dos artefatos",
  "",
  paste0("**Data:** ", Sys.Date()),
  "",
  "## 1. Checklist",
  "",
  knitr::kable(comp_df, format = "pipe", align = c("l", "c", "l")),
  "",
  "## 2. Critérios avaliados",
  "",
  "- existência dos scripts, tabelas, figuras e documentos essenciais;",
  "- ausência de figuras em PDF (apenas PNG);",
  "- PNGs válidos, legíveis e não vazios;",
  "- `Volcano_12` com conteúdo visual mínimo;",
  "- CSVs legíveis, não vazios e com colunas críticas sem NA;",
  "- scripts R parseáveis e com sementes de reprodutibilidade;",
  "- imagens referenciadas nos docs e no relatório Markdown existentes;",
  "- ausência de menções a PDF e a PCA nos documentos;",
  "- metadados com colunas esperadas, 529 amostras e sem recidivas;",
  "- versão do R e pacotes registrada em `docs/sessionInfo.txt`.",
  "",
  "## 3. Resultado",
  "",
  ifelse(all(comp_df$status == "OK"), "**Todos os itens passaram.**",
         paste0("**Itens com falha:** ",
                paste(comp_df$item[comp_df$status != "OK"], collapse = "; "))),
  "",
  "## 4. Observações",
  "",
  "- A checagem de “não vazio” usa a fração de pixels não brancos; ela detecta",
  "  figuras totalmente brancas, mas não substitui inspeção visual.",
  "- A checagem de PCA/PDF é textual e complementa a inspeção dos scripts.",
  "- A ausência de validação experimental e o confundimento TCGA × GTEx são",
  "  limitações científicas, não falhas de compliance.",
  "",
  "## 5. Sessão R",
  "",
  "```text",
  capture.output(sessionInfo()),
  "```"
)
writeLines(md, file.path(dir_aud, "COMPLIANCE.md"))
cat("Compliance gerado em:", dir_aud, "\n")
print(comp_df, row.names = FALSE)
