# =============================================================================
# etapa_08_analises_grupos_separados.R
# -----------------------------------------------------------------------------
# Análises separadas dos grupos alternativos do TCGA-LIHC:
#   1. Tumores recidivados (n = 2), incluindo expressão no painel LA;
#   2. Subtipos histológicos (HCC clássico, misto e fibrolamelar);
#   3. Grupos virais (HBV, HCV, HBV+HCV e Unknown);
#   4. Centros TCGA (TSS), incluindo sinal de lote no UMAP.
#
# Princípio: n pequeno e missingness alta => análises descritivas/exploratórias,
# sem inferência forte. As recidivas não entram na análise principal.
#
# Saídas: tabelas 23–35, figuras 28–35 e resumo em resultados/controle_qualidade/.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
dir_proc <- file.path(ROOT, "dados", "processed")
dir_tab  <- file.path(ROOT, "resultados", "tabelas")
dir_fig  <- file.path(ROOT, "resultados", "figuras")
dir_qc   <- file.path(ROOT, "resultados", "controle_qualidade")
for (d in c(dir_tab, dir_fig, dir_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

# -----------------------------------------------------------------------------
# 0. Dados
# -----------------------------------------------------------------------------
meta_all <- fread(file.path(dir_proc, "metadata_amostras.tsv"),
                  sep = "\t", data.table = FALSE, na.strings = c("", "NA"))
meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
              sep = "\t", data.table = FALSE, na.strings = c("", "NA"))

prim <- meta_all[meta_all$group == "TCGA_Tumor", ]
rec  <- meta_all[meta_all$group == "TCGA_Recurrent", ]
tcga <- meta[meta$group == "TCGA_Tumor", ]

tcga$sex <- factor(tcga$sex, levels = c("female", "male"))
tcga$stage_group <- factor(tcga$stage_group, levels = c("I-II", "III-IV"))
tcga$grade_group <- factor(tcga$grade_group, levels = c("G1-G2", "G3-G4"))
tcga$viral_status <- factor(tcga$viral_status,
                            levels = c("HBV", "HCV", "HBV+HCV", "Unknown"))
tcga$viral_known <- factor(
  ifelse(tcga$viral_status == "Unknown", NA, as.character(tcga$viral_status)),
  levels = c("HBV", "HCV", "HBV+HCV")
)
tcga$viral_flag <- factor(
  ifelse(tcga$viral_status == "Unknown", "Unknown", "Known"),
  levels = c("Known", "Unknown")
)
tcga$histology_label <- factor(
  c("Hepatocellular Carcinoma" = "HCC clássico",
    "Hepatocholangiocarcinoma (Mixed)" = "Misto",
    "Fibrolamellar Carcinoma" = "Fibrolamelar")[as.character(tcga$histology)],
  levels = c("HCC clássico", "Misto", "Fibrolamelar")
)
tcga$tss <- as.character(tcga$tss)

marcadores <- c("GPC3", "AFP", "KRT19", "EPCAM", "BIRC5",
                "VIM", "SPP1", "MMP9", "COL1A1")

carregar_expressao <- function(arquivo) {
  dt <- fread(arquivo, sep = "\t", data.table = FALSE)
  rownames(dt) <- dt$gene
  m <- as.matrix(dt[, setdiff(names(dt), "gene")])
  storage.mode(m) <- "numeric"
  m
}

expr_full <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))
expr_la   <- carregar_expressao(file.path(dir_proc, "expressao_LA_pathway.tsv.gz"))

# -----------------------------------------------------------------------------
# Utilidades estatísticas e de formatação
# -----------------------------------------------------------------------------
med <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
q1  <- function(x) if (sum(!is.na(x)) < 2) NA_real_ else as.numeric(quantile(x, 0.25, na.rm = TRUE))
q3  <- function(x) if (sum(!is.na(x)) < 2) NA_real_ else as.numeric(quantile(x, 0.75, na.rm = TRUE))
pct_true <- function(x) if (all(is.na(x))) NA_real_ else 100 * mean(x, na.rm = TRUE)
fmt_num <- function(x, digits = 1) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits, decimal.mark = ","))
fmt_pct <- function(n, total) {
  if (is.na(n) || total == 0) return("NA")
  sub("\\.", ",", sprintf("%d (%.1f%%)", n, 100 * n / total))
}
fmt_med_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  if (length(x) == 1) return(fmt_num(x, digits))
  paste0(fmt_num(med(x), digits), " (", fmt_num(q1(x), digits), "–", fmt_num(q3(x), digits), ")")
}

kw_p <- function(value, group) {
  d <- data.frame(value = value, group = group)
  d <- d[is.finite(d$value) & !is.na(d$group), ]
  d$group <- droplevels(factor(d$group))
  if (nrow(d) < 5 || nlevels(d$group) < 2) return(NA_real_)
  tryCatch(kruskal.test(value ~ group, data = d)$p.value, error = function(e) NA_real_)
}

cat_p <- function(x, y) {
  d <- data.frame(x = x, y = y)
  d <- d[!is.na(d$x) & !is.na(d$y), ]
  d$x <- droplevels(factor(d$x)); d$y <- droplevels(factor(d$y))
  tab <- table(d$x, d$y)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) {
    suppressWarnings(fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value)
  } else {
    suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
  }
}

summ_clin <- function(d, gvar) {
  d$grupo <- d[[gvar]]
  d %>%
    filter(!is.na(grupo)) %>%
    group_by(grupo) %>%
    summarise(
      n = n(),
      idade_mediana = med(age),
      idade_q1 = q1(age),
      idade_q3 = q3(age),
      masculino_pct = pct_true(sex == "male"),
      estadio_III_IV_pct = pct_true(stage_group == "III-IV"),
      grau_G3_G4_pct = pct_true(grade_group == "G3-G4"),
      afp_mediana = med(afp),
      afp_q1 = q1(afp),
      afp_q3 = q3(afp),
      fibrose_mediana = med(fibrosis_ishak),
      obitos = sum(OS == 1, na.rm = TRUE),
      os_mediana = med(OS_time),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
}

# =============================================================================
# 1. TUMORES RECIDIVADOS (n = 2)
# =============================================================================
rec_cols <- c("sample", "subject", "sample_type", "histology", "age", "sex",
              "stage", "grade", "afp", "fibrosis_ishak", "viral_status",
              "OS", "OS_time", "vital_status", "tss")
rec_out <- rec[, rec_cols]
write.csv(rec_out, file.path(dir_tab, "tabela_23_recidivas_clinico.csv"), row.names = FALSE)

prim_n <- nrow(prim)
comp <- data.frame(
  variavel = c("n", "Idade (anos)", "Sexo masculino", "Estágio I-II",
               "Grau G3-G4", "AFP (ng/mL), mediana (IQR)",
               "Fibrose Ishak, mediana", "Vírus detectável", "Óbito"),
  recidiva_1 = c(
    "1", fmt_num(rec$age[1]),
    ifelse(rec$sex[1] == "male", "sim", "não"),
    ifelse(rec$stage[1] %in% c("I", "II"), "sim", "não"),
    ifelse(rec$grade[1] %in% c("G3", "G4"), "sim", "não"),
    fmt_med_iqr(rec$afp[1]),
    fmt_num(rec$fibrosis_ishak[1]),
    ifelse(rec$viral_status[1] %in% c("HBV", "HCV", "HBV+HCV"),
           paste0("sim (", rec$viral_status[1], ")"), paste0("não (", rec$viral_status[1], ")")),
    ifelse(rec$OS[1] == 1, "sim", "não")
  ),
  recidiva_2 = c(
    "1", fmt_num(rec$age[2]),
    ifelse(rec$sex[2] == "male", "sim", "não"),
    ifelse(rec$stage[2] %in% c("I", "II"), "sim", "não"),
    ifelse(rec$grade[2] %in% c("G3", "G4"), "sim", "não"),
    fmt_med_iqr(rec$afp[2]),
    fmt_num(rec$fibrosis_ishak[2]),
    ifelse(rec$viral_status[2] %in% c("HBV", "HCV", "HBV+HCV"),
           paste0("sim (", rec$viral_status[2], ")"), paste0("não (", rec$viral_status[2], ")")),
    ifelse(rec$OS[2] == 1, "sim", "não")
  ),
  primarios = c(
    as.character(prim_n),
    fmt_med_iqr(prim$age),
    fmt_pct(sum(prim$sex == "male", na.rm = TRUE), prim_n),
    fmt_pct(sum(prim$stage %in% c("I", "II"), na.rm = TRUE), sum(!is.na(prim$stage))),
    fmt_pct(sum(prim$grade %in% c("G3", "G4"), na.rm = TRUE), sum(!is.na(prim$grade))),
    fmt_med_iqr(prim$afp),
    fmt_num(med(prim$fibrosis_ishak)),
    fmt_pct(sum(prim$viral_status %in% c("HBV", "HCV", "HBV+HCV"), na.rm = TRUE),
            prim_n),
    fmt_pct(sum(prim$OS == 1, na.rm = TRUE), sum(!is.na(prim$OS)))
  ),
  stringsAsFactors = FALSE
)
write.csv(comp, file.path(dir_tab, "tabela_24_recidivas_comparacao.csv"), row.names = FALSE)

# Expressão das recidivas no painel LA (a matriz completa não inclui -02)
samples_prim <- intersect(prim$sample, colnames(expr_la))
samples_rec  <- intersect(rec$sample, colnames(expr_la))
X <- expr_la[, c(samples_prim, samples_rec), drop = FALSE]
mu <- rowMeans(X[, samples_prim, drop = FALSE], na.rm = TRUE)
sdv <- apply(X[, samples_prim, drop = FALSE], 1, sd, na.rm = TRUE)
keep <- is.finite(mu) & is.finite(sdv) & sdv > 0
Xok <- X[keep, , drop = FALSE]
mu <- mu[keep]; sdv <- sdv[keep]
z <- sweep(sweep(Xok[, samples_rec, drop = FALSE], 1, mu), 1, sdv, "/")
pct <- sapply(samples_rec, function(s) {
  v <- Xok[, s]
  sapply(seq_along(v), function(i) {
    vals <- Xok[i, samples_prim]
    100 * mean(vals <= v[i], na.rm = TRUE)
  })
})
pct <- as.matrix(pct)
rownames(pct) <- rownames(Xok)

maxz <- apply(abs(z), 1, max, na.rm = TRUE)
top_genes <- names(sort(maxz, decreasing = TRUE))[1:10]
rec_expr <- do.call(rbind, lapply(samples_rec, function(s) {
  data.frame(
    sample = s, gene = top_genes,
    expressao = as.numeric(Xok[top_genes, s]),
    z = as.numeric(z[top_genes, s]),
    percentil_entre_primarios = as.numeric(pct[top_genes, s]),
    stringsAsFactors = FALSE
  )
}))
write.csv(rec_expr, file.path(dir_tab, "tabela_25_recidivas_expressao.csv"), row.names = FALSE)

prim_long <- data.frame(
  gene = rep(top_genes, each = length(samples_prim)),
  group = "Primário",
  expressao = as.numeric(Xok[top_genes, samples_prim])
)
rec_long <- data.frame(
  gene = rep(top_genes, times = length(samples_rec)),
  group = rep(ifelse(seq_along(samples_rec) == 1, "Recidiva 1", "Recidiva 2"),
              each = length(top_genes)),
  expressao = as.numeric(Xok[top_genes, samples_rec])
)
prim_long$gene <- factor(prim_long$gene, levels = top_genes)
rec_long$gene <- factor(rec_long$gene, levels = top_genes)
p_rec <- ggplot() +
  geom_violin(data = prim_long, aes(x = group, y = expressao),
              fill = "grey78", colour = NA, alpha = 0.65, trim = FALSE,
              scale = "width") +
  geom_boxplot(data = prim_long, aes(x = group, y = expressao),
               width = 0.12, outlier.shape = NA, alpha = 0.9,
               colour = "grey25") +
  geom_point(data = rec_long, aes(x = group, y = expressao, colour = group),
             size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_x_discrete(limits = c("Primário", "Recidiva 1", "Recidiva 2")) +
  scale_colour_manual(values = c("Recidiva 1" = "#C0392B",
                                 "Recidiva 2" = "#7D3C98")) +
  labs(title = "Expressão das recidivas no painel LA",
       subtitle = "Top 10 genes com maior |z| nas recidivas; violino = 369 tumores primários",
       x = NULL, y = "Expressão (log2 TPM + 0,001)", colour = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(dir_fig, "figura_28_recidivas_expressao.png"), p_rec,
       width = 11, height = 6.5, dpi = 300)

# Correlação global (painel LA) das recidivas com cada grupo
grupos_cor <- list(
  "Primário" = intersect(prim$sample, colnames(expr_la)),
  "Adjacente" = intersect(meta_all$sample[meta_all$group == "TCGA_Adjacent"], colnames(expr_la)),
  "Normal" = intersect(meta_all$sample[meta_all$group == "GTEx_Normal"], colnames(expr_la))
)
cor_summary <- list()
cor_long <- list()
for (s in samples_rec) {
  for (g in names(grupos_cor)) {
    ss <- grupos_cor[[g]]
    cors <- apply(expr_la[, ss, drop = FALSE], 2, function(v) {
      cor(expr_la[, s], v, method = "spearman", use = "pairwise.complete.obs")
    })
    cor_summary[[length(cor_summary) + 1]] <- data.frame(
      sample = s, grupo = g, n = length(cors),
      rho_mediana = median(cors, na.rm = TRUE),
      rho_maximo = max(cors, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    cor_long[[length(cor_long) + 1]] <- data.frame(
      sample = s, grupo = g, rho = as.numeric(cors)
    )
  }
}
cor_summary <- do.call(rbind, cor_summary)
cor_long <- do.call(rbind, cor_long)
write.csv(cor_summary, file.path(dir_tab, "tabela_26_recidivas_correlacao.csv"), row.names = FALSE)

cor_long$grupo <- factor(cor_long$grupo, levels = c("Primário", "Adjacente", "Normal"))
cor_long$sample_curto <- ifelse(cor_long$sample == samples_rec[1], "Recidiva 1", "Recidiva 2")
p_cor <- ggplot(cor_long, aes(x = grupo, y = rho, fill = grupo)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9, colour = "grey25") +
  facet_wrap(~ sample_curto) +
  scale_fill_manual(values = c("Primário" = "#C0392B", "Adjacente" = "#E67E22",
                               "Normal" = "#2471A3")) +
  labs(title = "Similaridade global das recidivas com os grupos",
       subtitle = "Correlação de Spearman por amostra, usando os 213 genes do painel LA",
       x = NULL, y = "Correlação de Spearman (ρ)", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")
ggsave(file.path(dir_fig, "figura_29_recidivas_correlacao.png"), p_cor,
       width = 8, height = 4.5, dpi = 300)

# =============================================================================
# 2. SUBTIPOS HISTOLÓGICOS
# =============================================================================
hist_clin <- summ_clin(tcga, "histology_label")
write.csv(hist_clin, file.path(dir_tab, "tabela_27_histologia_clinico.csv"), row.names = FALSE)

hist_tests <- list()
for (g in marcadores) {
  d <- data.frame(value = as.numeric(expr_full[g, tcga$sample]),
                  hist = tcga$histology_label)
  d <- d[is.finite(d$value) & !is.na(d$hist), ]
  p_kw <- kw_p(d$value, d$hist)
  d_mix <- d[d$hist %in% c("HCC clássico", "Misto"), ]
  d_mix$hist <- droplevels(d_mix$hist)
  p_mix <- if (nlevels(d_mix$hist) == 2 && min(table(d_mix$hist)) >= 2) {
    suppressWarnings(wilcox.test(value ~ hist, data = d_mix)$p.value)
  } else NA_real_
  d_fib <- d[d$hist %in% c("HCC clássico", "Fibrolamelar"), ]
  d_fib$hist <- droplevels(d_fib$hist)
  p_fib <- if (nlevels(d_fib$hist) == 2 && min(table(d_fib$hist)) >= 2) {
    suppressWarnings(wilcox.test(value ~ hist, data = d_fib)$p.value)
  } else NA_real_
  hist_tests[[length(hist_tests) + 1]] <- data.frame(
    gene = g, n = nrow(d),
    p_kruskal = p_kw,
    p_misto_vs_HCC = p_mix,
    p_fibrolamelar_vs_HCC = p_fib,
    stringsAsFactors = FALSE
  )
}
hist_tests <- do.call(rbind, hist_tests)
hist_tests$p_kruskal_BH <- p.adjust(hist_tests$p_kruskal, method = "BH")
hist_tests$p_misto_BH <- p.adjust(hist_tests$p_misto_vs_HCC, method = "BH")
hist_tests$p_fibro_BH <- p.adjust(hist_tests$p_fibrolamelar_vs_HCC, method = "BH")
write.csv(hist_tests, file.path(dir_tab, "tabela_28_histologia_testes.csv"), row.names = FALSE)

p_hist_age <- ggplot(tcga[!is.na(tcga$age) & !is.na(tcga$histology_label), ],
                     aes(x = histology_label, y = age, fill = histology_label)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.35) +
  scale_fill_manual(values = c("HCC clássico" = "#C0392B", "Misto" = "#7D3C98",
                               "Fibrolamelar" = "#F1C40F")) +
  labs(title = "Idade por subtipo histológico", x = NULL, y = "Idade (anos)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

p_hist_afp <- ggplot(tcga[!is.na(tcga$afp_log2) & !is.na(tcga$histology_label), ],
                     aes(x = histology_label, y = afp_log2, fill = histology_label)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.35) +
  scale_fill_manual(values = c("HCC clássico" = "#C0392B", "Misto" = "#7D3C98",
                               "Fibrolamelar" = "#F1C40F")) +
  labs(title = "AFP por subtipo histológico", x = NULL,
       y = "log2(AFP + 1)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

p_hist <- p_hist_age | p_hist_afp
ggsave(file.path(dir_fig, "figura_30_histologia_clinico.png"), p_hist,
       width = 10, height = 4.8, dpi = 300)

# =============================================================================
# 3. GRUPOS VIRAIS
# =============================================================================
viral_clin <- summ_clin(tcga, "viral_status")
write.csv(viral_clin, file.path(dir_tab, "tabela_29_viral_clinico.csv"), row.names = FALSE)

viral_tests <- list()
add_viral <- function(comparacao, variavel, teste, n, p) {
  viral_tests[[length(viral_tests) + 1]] <<- data.frame(
    comparacao = comparacao, variavel = variavel, teste = teste,
    n = n, p_valor = p, stringsAsFactors = FALSE
  )
}

# Entre os grupos virais conhecidos (Unknown excluído)
for (v in c("age", "afp_log2", "fibrosis_ishak")) {
  d <- tcga[!is.na(tcga$viral_known) & !is.na(tcga[[v]]), ]
  if (nrow(d) >= 10 && nlevels(droplevels(d$viral_known)) >= 2) {
    p <- suppressWarnings(kruskal.test(d[[v]] ~ droplevels(d$viral_known))$p.value)
    add_viral("HBV vs HCV vs HBV+HCV", v, "Kruskal-Wallis", nrow(d), p)
  }
}
for (v in c("sex", "stage_group", "grade_group")) {
  d <- tcga[!is.na(tcga$viral_known) & !is.na(tcga[[v]]), ]
  if (nrow(d) >= 10) {
    p <- cat_p(d$viral_known, d[[v]])
    add_viral("HBV vs HCV vs HBV+HCV", v, "Qui-quadrado/Fisher", nrow(d), p)
  }
}

# Conhecido vs Unknown (avalia se a missingness é informativa)
for (v in c("age", "afp_log2", "fibrosis_ishak")) {
  d <- tcga[!is.na(tcga[[v]]), ]
  if (nrow(d) >= 10 && nlevels(droplevels(d$viral_flag)) == 2) {
    p <- suppressWarnings(wilcox.test(d[[v]] ~ droplevels(d$viral_flag))$p.value)
    add_viral("Conhecido vs Unknown", v, "Wilcoxon", nrow(d), p)
  }
}
for (v in c("sex", "stage_group", "grade_group")) {
  d <- tcga[!is.na(tcga[[v]]), ]
  if (nrow(d) >= 10) {
    p <- cat_p(d$viral_flag, d[[v]])
    add_viral("Conhecido vs Unknown", v, "Qui-quadrado/Fisher", nrow(d), p)
  }
}
viral_tests <- do.call(rbind, viral_tests)
viral_tests$p_ajustado_BH <- p.adjust(viral_tests$p_valor, method = "BH")
viral_tests$significativo_BH <- ifelse(viral_tests$p_ajustado_BH < 0.05, "Sim", "Não")
write.csv(viral_tests, file.path(dir_tab, "tabela_30_viral_testes.csv"), row.names = FALSE)

viral_markers <- list()
viral_p <- list()
for (g in marcadores) {
  for (v in levels(tcga$viral_status)) {
    vals <- as.numeric(expr_full[g, tcga$sample[tcga$viral_status == v]])
    vals <- vals[is.finite(vals)]
    viral_markers[[length(viral_markers) + 1]] <- data.frame(
      gene = g, viral_status = v, n = length(vals),
      mediana = med(vals), q1 = q1(vals), q3 = q3(vals),
      stringsAsFactors = FALSE
    )
  }
  p <- kw_p(as.numeric(expr_full[g, tcga$sample]), tcga$viral_known)
  viral_p[[length(viral_p) + 1]] <- data.frame(gene = g, p_kw_conhecidos = p)
}
viral_markers <- do.call(rbind, viral_markers)
viral_p <- do.call(rbind, viral_p)
viral_p$p_BH <- p.adjust(viral_p$p_kw_conhecidos, method = "BH")
viral_markers <- merge(viral_markers, viral_p, by = "gene", all.x = TRUE)
viral_markers$viral_status <- factor(viral_markers$viral_status,
                                     levels = levels(tcga$viral_status))
viral_markers <- viral_markers[order(viral_markers$gene, viral_markers$viral_status), ]
write.csv(viral_markers, file.path(dir_tab, "tabela_31_viral_marcadores.csv"), row.names = FALSE)

p_viral_afp <- ggplot(tcga[!is.na(tcga$afp_log2) & !is.na(tcga$viral_status), ],
                      aes(x = viral_status, y = afp_log2, fill = viral_status)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.7, alpha = 0.3) +
  scale_fill_manual(values = c("HBV" = "#C0392B", "HCV" = "#E67E22",
                               "HBV+HCV" = "#7D3C98", "Unknown" = "#95A5A6")) +
  labs(title = "AFP por grupo viral",
       subtitle = "Unknown mostrado para transparência; testes entre os grupos conhecidos",
       x = NULL, y = "log2(AFP + 1)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(dir_fig, "figura_31_viral_clinico.png"), p_viral_afp,
       width = 7, height = 4.5, dpi = 300)

viral_long <- as.data.frame(t(expr_full[marcadores, tcga$sample, drop = FALSE]))
viral_long$sample <- rownames(viral_long)
viral_long <- merge(viral_long, tcga[, c("sample", "viral_status")], by = "sample")
viral_long <- pivot_longer(viral_long, cols = all_of(marcadores),
                           names_to = "gene", values_to = "expressao")
viral_long$gene <- factor(viral_long$gene, levels = marcadores)
p_viral_markers <- ggplot(viral_long, aes(x = viral_status, y = expressao,
                                          fill = viral_status)) +
  geom_violin(trim = FALSE, alpha = 0.5, colour = NA, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9,
               colour = "grey25") +
  geom_jitter(width = 0.07, size = 0.45, alpha = 0.25) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("HBV" = "#C0392B", "HCV" = "#E67E22",
                               "HBV+HCV" = "#7D3C98", "Unknown" = "#95A5A6")) +
  labs(title = "Marcadores por grupo viral — TCGA tumor",
       subtitle = "Testes de Kruskal-Wallis apenas entre HBV, HCV e HBV+HCV",
       x = NULL, y = "Expressão (log2 TPM + 0,001)", fill = "Vírus") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(dir_fig, "figura_32_viral_marcadores.png"), p_viral_markers,
       width = 11, height = 8, dpi = 300)

# =============================================================================
# 4. CENTROS TCGA (TSS)
# =============================================================================
tss_summary <- summ_clin(tcga, "tss")
write.csv(tss_summary, file.path(dir_tab, "tabela_32_tss_resumo.csv"), row.names = FALSE)

top_tss <- tss_summary$grupo[tss_summary$n >= 10]
top_tss <- as.character(top_tss)
tcga$tss_group <- factor(ifelse(tcga$tss %in% top_tss, tcga$tss, "Outros"),
                         levels = c(top_tss, "Outros"))

tss_tests <- list()
add_tss <- function(variavel, teste, n, p) {
  tss_tests[[length(tss_tests) + 1]] <<- data.frame(
    variavel = variavel, teste = teste, n = n, p_valor = p,
    stringsAsFactors = FALSE
  )
}
d_top <- tcga[tcga$tss %in% top_tss, ]
for (v in c("age", "afp_log2", "fibrosis_ishak")) {
  d <- d_top[!is.na(d_top[[v]]), ]
  add_tss(v, "Kruskal-Wallis", nrow(d), kw_p(d[[v]], droplevels(factor(d$tss))))
}
for (v in c("sex", "stage_group", "grade_group")) {
  d <- d_top[!is.na(d_top[[v]]), ]
  add_tss(v, "Qui-quadrado/Fisher", nrow(d), cat_p(d$tss, d[[v]]))
}
for (g in marcadores) {
  d <- d_top
  vals <- as.numeric(expr_full[g, d$sample])
  add_tss(paste0("expressao_", g), "Kruskal-Wallis", length(vals),
          kw_p(vals, droplevels(factor(d$tss))))
}
tss_tests <- do.call(rbind, tss_tests)
tss_tests$p_ajustado_BH <- p.adjust(tss_tests$p_valor, method = "BH")
write.csv(tss_tests, file.path(dir_tab, "tabela_33_tss_testes.csv"), row.names = FALSE)

expr_top <- as.data.frame(t(expr_full[marcadores, d_top$sample, drop = FALSE]))
expr_top$sample <- rownames(expr_top)
d_top_expr <- merge(d_top[, c("sample", "tss")], expr_top, by = "sample")
tss_marker <- d_top_expr %>%
  group_by(tss) %>%
  summarise(n = n(), across(all_of(marcadores), med), .groups = "drop") %>%
  arrange(desc(n))
write.csv(tss_marker, file.path(dir_tab, "tabela_34_tss_marcadores.csv"), row.names = FALSE)

# UMAP: sinal de lote associado ao TSS
umap_t <- fread(file.path(dir_tab, "tabela_09_UMAP_coordenadas.csv"),
                data.table = FALSE, na.strings = c("", "NA"))
umap_t <- umap_t[umap_t$group == "TCGA_Tumor", ]
umap_t$tss_group <- factor(ifelse(umap_t$tss %in% top_tss, umap_t$tss, "Outros"),
                           levels = c(top_tss, "Outros"))
umap_tests <- list()
for (pc in c("UMAP1", "UMAP2")) {
  d <- umap_t[!is.na(umap_t[[pc]]) & !is.na(umap_t$tss_group), ]
  fit <- aov(d[[pc]] ~ d$tss_group)
  a <- summary(fit)[[1]]
  r2 <- a[["Sum Sq"]][1] / sum(a[["Sum Sq"]])
  p <- a[["Pr(>F)"]][1]
  umap_tests[[length(umap_tests) + 1]] <- data.frame(
    coordenada = pc, n = nrow(d),
    r2_tss = r2, p_aov = p,
    p_kruskal = suppressWarnings(kruskal.test(d[[pc]] ~ d$tss_group)$p.value),
    stringsAsFactors = FALSE
  )
}
umap_tests <- do.call(rbind, umap_tests)
write.csv(umap_tests, file.path(dir_tab, "tabela_35_tss_umap.csv"), row.names = FALSE)

# Figura: distribuição de amostras por TSS (top 14 + Outros)
tss_plot <- tss_summary %>%
  mutate(tss2 = ifelse(row_number() <= 14, grupo, "Outros")) %>%
  group_by(tss2) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(desc(n)) %>%
  mutate(tss2 = factor(tss2, levels = tss2))
p_tss_dist <- ggplot(tss_plot, aes(x = tss2, y = n,
                                    fill = tss2 %in% top_tss)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#95A5A6")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Amostras por centro TCGA (TSS)",
       subtitle = paste0("38 centros; ", length(top_tss),
                         " com n >= 10. Vermelho = centros incluídos nos testes."),
       x = NULL, y = "Número de tumores") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(dir_fig, "figura_33_tss_distribuicao.png"), p_tss_dist,
       width = 9, height = 5, dpi = 300)

# Figura: UMAP colorido por TSS (top centros + Outros)
cores_tss <- c(setNames(brewer.pal(max(3, length(top_tss)), "Dark2"), top_tss),
               "Outros" = "grey75")
p_tss_umap <- ggplot(umap_t, aes(x = UMAP1, y = UMAP2, colour = tss_group)) +
  geom_point(size = 1.1, alpha = 0.7) +
  scale_colour_manual(values = cores_tss) +
  labs(title = "UMAP dos tumores colorido por centro TCGA (TSS)",
       subtitle = "Centros com n >= 10 em cores; demais agregados em 'Outros'",
       x = "UMAP 1", y = "UMAP 2", colour = "TSS") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")
ggsave(file.path(dir_fig, "figura_34_tss_umap.png"), p_tss_umap,
       width = 8, height = 6, dpi = 300)

# Heatmap: medianas de marcadores por TSS (top centros)
mat <- as.matrix(tss_marker[, marcadores])
rownames(mat) <- tss_marker$tss
mat_z <- t(scale(mat))  # genes nas linhas, centros nas colunas
mat_z[is.na(mat_z)] <- 0
ann_col <- data.frame(n = tss_marker$n, row.names = tss_marker$tss)
png(file.path(dir_fig, "figura_35_tss_marcadores_heatmap.png"),
    width = 10, height = 5.5, units = "in", res = 300)
pheatmap(mat_z,
         annotation_col = ann_col,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-2, 2, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         cellwidth = 30, cellheight = 20,
         fontsize = 9, fontsize_row = 9, fontsize_col = 9,
         border_color = "grey90",
         main = "Medianas de marcadores por TSS (z-score entre centros)")
dev.off()

# =============================================================================
# 5. Resumo
# =============================================================================
resumo <- c(
  "GRUPOS SEPARADOS — RESUMO",
  "=========================",
  paste("Tumores primários:", nrow(tcga), "| recidivas:", nrow(rec)),
  "",
  "Recidivas — clínica:",
  paste(capture.output(print(rec_out, row.names = FALSE)), collapse = "\n"),
  "",
  "Recidivas — genes mais extremos no painel LA:",
  paste(capture.output(print(rec_expr, row.names = FALSE)), collapse = "\n"),
  "",
  "Recidivas — correlação global:",
  paste(capture.output(print(cor_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "Histologia — clínica:",
  paste(capture.output(print(hist_clin, row.names = FALSE)), collapse = "\n"),
  "",
  "Histologia — testes de marcadores:",
  paste(capture.output(print(hist_tests, row.names = FALSE)), collapse = "\n"),
  "",
  "Vírus — clínica:",
  paste(capture.output(print(viral_clin, row.names = FALSE)), collapse = "\n"),
  "",
  "Vírus — testes:",
  paste(capture.output(print(viral_tests, row.names = FALSE)), collapse = "\n"),
  "",
  "TSS — resumo (top 15):",
  paste(capture.output(print(head(tss_summary, 15), row.names = FALSE)), collapse = "\n"),
  "",
  "TSS — testes:",
  paste(capture.output(print(tss_tests, row.names = FALSE)), collapse = "\n"),
  "",
  "TSS — UMAP:",
  paste(capture.output(print(umap_tests, row.names = FALSE)), collapse = "\n"),
  "",
  "Sessão R:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(resumo, file.path(dir_qc, "resumo_grupos_separados.txt"))
cat(paste(resumo, collapse = "\n"), "\n")
