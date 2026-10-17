# =============================================================================
# 05_fenotipos_avancados.R
# -----------------------------------------------------------------------------
# Exploração ampliada dos fenótipos clínicos do TCGA-LIHC:
#   - missingness e qualidade das variáveis clínicas;
#   - associações entre sexo, idade, estágio, grau, AFP, fibrose, vírus,
#     Child-Pugh, histologia, invasão vascular e raça/etnia;
#   - sobrevida global (Kaplan-Meier e log-rank);
#   - modelos de Cox univariados e multivariados;
#   - verificação da suposição de riscos proporcionais (cox.zph);
#   - correlações entre variáveis clínicas.
#
# Saídas: tabelas 13–20, figuras 15–25 e resumo em resultados/qc/.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(broom)
  library(patchwork)
  library(pheatmap)
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
dir_qc   <- file.path(ROOT, "resultados", "qc")
for (d in c(dir_tab, dir_fig, dir_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

# -----------------------------------------------------------------------------
# 0. Dados e recodificações
# -----------------------------------------------------------------------------
meta <- fread(file.path(dir_proc, "metadata_analise.tsv"),
              sep = "\t", data.table = FALSE, na.strings = c("", "NA"))
tcga <- meta[meta$group == "TCGA_Tumor", ]

tcga$stage_group <- factor(tcga$stage_group, levels = c("I-II", "III-IV"))
tcga$grade_group <- factor(tcga$grade_group, levels = c("G1-G2", "G3-G4"))
tcga$sex <- factor(tcga$sex, levels = c("female", "male"))
tcga$age_group <- factor(tcga$age_group, levels = c("<65", ">=65"))
tcga$afp_high <- factor(tcga$afp_high, levels = c("<400", ">=400"))
tcga$fibrosis_advanced <- factor(tcga$fibrosis_advanced, levels = c("<5", ">=5"))
tcga$child_pugh_group <- factor(tcga$child_pugh_group, levels = c("A", "B-C"))
tcga$histology_group <- factor(tcga$histology_group, levels = c("HCC", "Other"))
tcga$vascular_invasion <- factor(tcga$vascular_invasion,
                                 levels = c("None/Unknown", "Micro", "Macro"))
tcga$viral_status <- factor(tcga$viral_status,
                            levels = c("None/Other", "HBV", "HCV", "HBV+HCV", "Unknown"))
tcga$race_group <- factor(tcga$race_group)
tcga$ethnicity_group <- factor(tcga$ethnicity_group)
# Para testes, excluir "Unknown" (não é uma categoria biológica)
tcga$viral_known <- factor(
  ifelse(tcga$viral_status == "Unknown", NA, as.character(tcga$viral_status)),
  levels = c("None/Other", "HBV", "HCV", "HBV+HCV")
)
tcga$race_known <- factor(ifelse(tcga$race_group == "Unknown", NA,
                                 as.character(tcga$race_group)))
tcga$ethnicity_known <- factor(ifelse(tcga$ethnicity_group == "Unknown", NA,
                                      as.character(tcga$ethnicity_group)))

# Variáveis numéricas derivadas
tcga$stage_num <- as.numeric(factor(tcga$stage, levels = c("I", "II", "III", "IV")))
tcga$grade_num <- as.numeric(factor(tcga$grade, levels = c("G1", "G2", "G3", "G4")))

# -----------------------------------------------------------------------------
# 1. Missingness clínico
# -----------------------------------------------------------------------------
vars_miss <- c("age", "sex", "stage", "grade", "afp", "afp_high",
               "fibrosis_ishak", "fibrosis_advanced", "viral_status",
               "child_pugh", "vascular_invasion", "histology",
               "race", "ethnicity", "OS", "OS_time")

calc_missing <- function(x) {
  if (is.factor(x) || is.character(x)) {
    is.na(x) | x == "" | x == "Unknown" | x == "not reported"
  } else {
    is.na(x)
  }
}
miss <- data.frame(
  variavel = vars_miss,
  n_na = sapply(tcga[vars_miss], function(x) sum(calc_missing(x))),
  n_total = nrow(tcga)
)
miss$percentual_na <- round(100 * miss$n_na / miss$n_total, 1)
miss <- miss[order(-miss$percentual_na), ]
write.csv(miss, file.path(dir_tab, "13_fenotipos_missingness.csv"), row.names = FALSE)

p_miss <- ggplot(miss, aes(x = reorder(variavel, percentual_na), y = percentual_na)) +
  geom_col(fill = "#2E86C1", width = 0.72) +
  geom_text(aes(label = paste0(percentual_na, "%")), hjust = -0.15, size = 3.2) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(title = "Missingness das variáveis clínicas — TCGA-LIHC (tumor)",
       subtitle = paste0("n = ", nrow(tcga), " tumores primários"),
       x = NULL, y = "% ausente/desconhecido") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(dir_fig, "16_missingness_clinico.png"), p_miss,
       width = 8, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# 2. Testes de associação
# -----------------------------------------------------------------------------
teste_cat <- function(d, v1, v2, lab1, lab2) {
  d <- d[!is.na(d[[v1]]) & !is.na(d[[v2]]), ]
  if (nrow(d) < 10) return(NULL)
  tab <- table(d[[v1]], d[[v2]])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  exp_tab <- suppressWarnings(chisq.test(tab)$expected)
  teste <- if (any(exp_tab < 5)) "Fisher (simulado)" else "Qui-quadrado"
  p <- if (teste == "Fisher (simulado)") {
    fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value
  } else {
    chisq.test(tab, correct = FALSE)$p.value
  }
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  V <- sqrt(as.numeric(chi) / (sum(tab) * (min(nrow(tab), ncol(tab)) - 1)))
  data.frame(
    variavel_1 = lab1, variavel_2 = lab2,
    n = sum(tab), teste = teste,
    estatistica = round(as.numeric(chi), 3),
    p_valor = p,
    cramers_v = round(V, 3),
    stringsAsFactors = FALSE
  )
}

teste_cont <- function(d, v, g, lab_v, lab_g) {
  d <- d[!is.na(d[[v]]) & !is.na(d[[g]]), ]
  d[[g]] <- droplevels(factor(d[[g]]))
  if (nrow(d) < 10 || nlevels(d[[g]]) < 2) return(NULL)
  if (nlevels(d[[g]]) == 2) {
    teste <- "Wilcoxon"
    p <- wilcox.test(d[[v]] ~ d[[g]])$p.value
    stat <- as.numeric(wilcox.test(d[[v]] ~ d[[g]])$statistic)
  } else {
    teste <- "Kruskal-Wallis"
    p <- kruskal.test(d[[v]] ~ d[[g]])$p.value
    stat <- as.numeric(kruskal.test(d[[v]] ~ d[[g]])$statistic)
  }
  data.frame(
    variavel_1 = lab_v, variavel_2 = lab_g,
    n = nrow(d), teste = teste,
    estatistica = round(stat, 3),
    p_valor = p,
    cramers_v = NA_real_,
    stringsAsFactors = FALSE
  )
}

testes <- list()
add <- function(x) if (!is.null(x)) testes[[length(testes) + 1]] <<- x

# Categóricas
add(teste_cat(tcga, "sex", "stage_group", "Sexo", "Estágio (I-II vs III-IV)"))
add(teste_cat(tcga, "sex", "grade_group", "Sexo", "Grau (G1-G2 vs G3-G4)"))
add(teste_cat(tcga, "afp_high", "stage_group", "AFP (>=400)", "Estágio (I-II vs III-IV)"))
add(teste_cat(tcga, "afp_high", "grade_group", "AFP (>=400)", "Grau (G1-G2 vs G3-G4)"))
add(teste_cat(tcga, "fibrosis_advanced", "stage_group", "Fibrose avançada (Ishak >=5)", "Estágio"))
add(teste_cat(tcga, "fibrosis_advanced", "grade_group", "Fibrose avançada (Ishak >=5)", "Grau"))
add(teste_cat(tcga, "viral_known", "stage_group", "Vírus (marcadores; Unknown excluído)", "Estágio"))
add(teste_cat(tcga, "viral_known", "grade_group", "Vírus (marcadores; Unknown excluído)", "Grau"))
add(teste_cat(tcga, "child_pugh_group", "stage_group", "Child-Pugh (A vs B-C)", "Estágio"))
add(teste_cat(tcga, "child_pugh_group", "grade_group", "Child-Pugh (A vs B-C)", "Grau"))
add(teste_cat(tcga, "vascular_invasion", "stage_group", "Invasão vascular", "Estágio"))
add(teste_cat(tcga, "vascular_invasion", "grade_group", "Invasão vascular", "Grau"))
add(teste_cat(tcga, "histology_group", "stage_group", "Histologia (HCC vs outro)", "Estágio"))
add(teste_cat(tcga, "race_known", "stage_group", "Raça (Unknown excluído)", "Estágio"))
add(teste_cat(tcga, "race_known", "grade_group", "Raça (Unknown excluído)", "Grau"))

# Contínuas
add(teste_cont(tcga, "age", "sex", "Idade", "Sexo"))
add(teste_cont(tcga, "age", "stage_group", "Idade", "Estágio"))
add(teste_cont(tcga, "age", "grade_group", "Idade", "Grau"))
add(teste_cont(tcga, "age", "viral_known", "Idade", "Vírus (Unknown excluído)"))
add(teste_cont(tcga, "afp_log2", "sex", "log2(AFP+1)", "Sexo"))
add(teste_cont(tcga, "afp_log2", "stage_group", "log2(AFP+1)", "Estágio"))
add(teste_cont(tcga, "afp_log2", "grade_group", "log2(AFP+1)", "Grau"))
add(teste_cont(tcga, "fibrosis_ishak", "stage_group", "Fibrose Ishak", "Estágio"))
add(teste_cont(tcga, "fibrosis_ishak", "grade_group", "Fibrose Ishak", "Grau"))

testes_df <- do.call(rbind, testes)
testes_df$p_ajustado_BH <- p.adjust(testes_df$p_valor, method = "BH")
testes_df$significativo <- ifelse(testes_df$p_ajustado_BH < 0.05, "Sim", "Não")
testes_df <- testes_df[order(testes_df$p_valor), ]
write.csv(testes_df, file.path(dir_tab, "14_fenotipos_testes.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. AFP por subgrupo
# -----------------------------------------------------------------------------
summ_cont <- function(d, v, g, lab_v, lab_g) {
  d <- d[!is.na(d[[v]]) & !is.na(d[[g]]), ]
  d %>% group_by(.data[[g]]) %>%
    summarise(
      variavel = lab_v,
      grupo = as.character(dplyr::cur_group()[[1]]),
      n = n(),
      mediana = median(.data[[v]], na.rm = TRUE),
      q1 = quantile(.data[[v]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[v]], 0.75, na.rm = TRUE),
      media = mean(.data[[v]], na.rm = TRUE),
      dp = sd(.data[[v]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(variavel, grupo, n, mediana, q1, q3, media, dp)
}
afp_tab <- rbind(
  summ_cont(tcga, "afp", "stage_group", "AFP", "Estágio"),
  summ_cont(tcga, "afp", "grade_group", "AFP", "Grau"),
  summ_cont(tcga, "afp", "sex", "AFP", "Sexo"),
  summ_cont(tcga, "afp", "viral_known", "AFP", "Vírus (Unknown excluído)")
)
write.csv(afp_tab, file.path(dir_tab, "15_AFP_por_subgrupo.csv"), row.names = FALSE)

p_afp1 <- ggplot(tcga[!is.na(tcga$afp_log2) & !is.na(tcga$stage_group), ],
                 aes(x = stage_group, y = afp_log2, fill = stage_group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5, alpha = 0.85) +
  geom_jitter(width = 0.14, size = 0.8, alpha = 0.35) +
  scale_fill_manual(values = c("I-II" = "#9BC2E6", "III-IV" = "#1F4E79")) +
  labs(title = "AFP por estágio", x = NULL, y = "log2(AFP + 1)", fill = "Estágio") +
  theme_bw(base_size = 11) + theme(legend.position = "none",
                                   panel.grid.minor = element_blank())

p_afp2 <- ggplot(tcga[!is.na(tcga$afp_log2) & !is.na(tcga$grade_group), ],
                 aes(x = grade_group, y = afp_log2, fill = grade_group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5, alpha = 0.85) +
  geom_jitter(width = 0.14, size = 0.8, alpha = 0.35) +
  scale_fill_manual(values = c("G1-G2" = "#F5B7B1", "G3-G4" = "#922B21")) +
  labs(title = "AFP por grau", x = NULL, y = "log2(AFP + 1)", fill = "Grau") +
  theme_bw(base_size = 11) + theme(legend.position = "none",
                                   panel.grid.minor = element_blank())

p_afp <- p_afp1 | p_afp2
ggsave(file.path(dir_fig, "15_AFP_por_estagio_grau.png"), p_afp,
       width = 9, height = 4.5, dpi = 300)

# -----------------------------------------------------------------------------
# 4. Sobrevida global (Kaplan-Meier)
# -----------------------------------------------------------------------------
surv_df <- tcga[!is.na(tcga$OS_time) & tcga$OS_time > 0 & !is.na(tcga$OS), ]
surv_df$OS <- as.integer(surv_df$OS)

logrank_tab <- list()
km_plot <- function(var, titulo, arquivo, palette = "Dark2") {
  d <- surv_df[!is.na(surv_df[[var]]), ]
  d[[var]] <- droplevels(factor(d[[var]]))
  if (nlevels(d[[var]]) < 2) return(NULL)
  fml <- as.formula(paste("Surv(OS_time, OS) ~", var))
  fit <- survfit(fml, data = d)
  fit$call$formula <- fml
  sd <- survdiff(fml, data = d)
  p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  pobj <- ggsurvplot(
    fit, data = d, pval = TRUE, pval.method = FALSE,
    conf.int = FALSE, risk.table = TRUE, risk.table.height = 0.25,
    palette = palette, xlab = "Tempo (dias)", ylab = "Sobrevida global",
    legend.title = "", ggtheme = theme_bw(base_size = 11),
    tables.theme = theme_cleantable()
  )
  png(file.path(dir_fig, paste0(arquivo, ".png")), width = 8, height = 6,
      units = "in", res = 300)
  print(pobj)
  dev.off()
  data.frame(
    variavel = var, n = nrow(d), eventos = sum(d$OS),
    logrank_chisq = round(sd$chisq, 3), logrank_p = p,
    stringsAsFactors = FALSE
  )
}

for (spec in list(
  list("stage_group", "Sobrevida global por estágio", "17_KM_OS_estagio"),
  list("grade_group", "Sobrevida global por grau", "18_KM_OS_grau"),
  list("sex", "Sobrevida global por sexo", "19_KM_OS_sexo"),
  list("age_group", "Sobrevida global por idade", "20_KM_OS_idade"),
  list("afp_high", "Sobrevida global por AFP", "21_KM_OS_AFP"),
  list("fibrosis_advanced", "Sobrevida global por fibrose", "22_KM_OS_fibrose")
)) {
  x <- km_plot(spec[[1]], spec[[2]], spec[[3]])
  if (!is.null(x)) logrank_tab[[length(logrank_tab) + 1]] <- x
}
logrank_df <- do.call(rbind, logrank_tab)
write.csv(logrank_df, file.path(dir_tab, "16_KM_logrank.csv"), row.names = FALSE)

# Painel sem tabelas de risco (4 curvas)
painel_vars <- list(
  c("stage_group", "Estágio"),
  c("grade_group", "Grau"),
  c("sex", "Sexo"),
  c("age_group", "Idade")
)
plots_painel <- lapply(painel_vars, function(v) {
  d <- surv_df[!is.na(surv_df[[v[1]]]), ]
  d[[v[1]]] <- droplevels(factor(d[[v[1]]]))
  fml <- as.formula(paste("Surv(OS_time, OS) ~", v[1]))
  fit <- survfit(fml, data = d)
  fit$call$formula <- fml
  ggsurvplot(fit, data = d, pval = TRUE, conf.int = FALSE,
             legend.title = "", xlab = "Tempo (dias)", ylab = "Sobrevida global",
             palette = "Dark2", ggtheme = theme_bw(base_size = 10))$plot +
    labs(title = v[2])
})
p_painel <- wrap_plots(plots_painel, ncol = 2)
ggsave(file.path(dir_fig, "23_KM_OS_painel.png"), p_painel,
       width = 11, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Cox univariado e multivariado
# -----------------------------------------------------------------------------
cox_univ <- function(var, label) {
  d <- surv_df[!is.na(surv_df[[var]]), ]
  if (nrow(d) < 20) return(NULL)
  fit <- coxph(as.formula(paste("Surv(OS_time, OS) ~", var)), data = d)
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  td$variavel <- label
  td$n <- nrow(d)
  td$eventos <- sum(d$OS)
  td[, c("variavel", "term", "estimate", "conf.low", "conf.high",
         "p.value", "n", "eventos")]
}
cox_univ_list <- list()
for (spec in list(
  c("stage_group", "Estágio"), c("grade_group", "Grau"),
  c("sex", "Sexo"), c("age", "Idade"), c("age_group", "Idade (categorizada)"),
  c("afp_high", "AFP >=400"), c("afp_log2", "log2(AFP+1)"),
  c("fibrosis_advanced", "Fibrose avançada"),
  c("viral_known", "Vírus (Unknown excluído)"), c("child_pugh_group", "Child-Pugh"),
  c("vascular_invasion", "Invasão vascular"),
  c("histology_group", "Histologia")
)) {
  x <- cox_univ(spec[1], spec[2])
  if (!is.null(x)) cox_univ_list[[length(cox_univ_list) + 1]] <- x
}
cox_univ_df <- do.call(rbind, cox_univ_list)
write.csv(cox_univ_df, file.path(dir_tab, "17_cox_univariado.csv"), row.names = FALSE)

# Modelo A: estágio + grau + sexo + idade (n ≈ 342)
varsA <- c("OS_time", "OS", "stage_group", "grade_group", "sex", "age")
dA <- surv_df[complete.cases(surv_df[, varsA]), ]
fitA <- coxph(Surv(OS_time, OS) ~ stage_group + grade_group + sex + age, data = dA)
tidyA <- broom::tidy(fitA, exponentiate = TRUE, conf.int = TRUE)
tidyA$modelo <- "Modelo A: estágio + grau + sexo + idade"
tidyA$n <- nrow(dA); tidyA$eventos <- sum(dA$OS)

# Modelo B: A + AFP + fibrose (n ≈ 180)
varsB <- c(varsA, "afp_high", "fibrosis_advanced")
dB <- surv_df[complete.cases(surv_df[, varsB]), ]
fitB <- coxph(Surv(OS_time, OS) ~ stage_group + grade_group + sex + age +
                afp_high + fibrosis_advanced, data = dB)
tidyB <- broom::tidy(fitB, exponentiate = TRUE, conf.int = TRUE)
tidyB$modelo <- "Modelo B: A + AFP + fibrose"
tidyB$n <- nrow(dB); tidyB$eventos <- sum(dB$OS)

cox_multi <- rbind(tidyA, tidyB)
write.csv(cox_multi, file.path(dir_tab, "18_cox_multivariado.csv"), row.names = FALSE)

# Riscos proporcionais
phA <- cox.zph(fitA)
phB <- cox.zph(fitB)
ph_df <- rbind(
  data.frame(modelo = "Modelo A", as.data.frame(phA$table), check.names = FALSE),
  data.frame(modelo = "Modelo B", as.data.frame(phB$table), check.names = FALSE)
)
ph_df$termo <- rownames(ph_df)
rownames(ph_df) <- NULL
ph_df <- ph_df[, c("modelo", "termo", "chisq", "df", "p")]
write.csv(ph_df, file.path(dir_tab, "20_cox_ph_assumption.csv"), row.names = FALSE)

# Forest plot
rotulos <- c(
  "stage_groupIII-IV" = "Estágio III-IV (vs I-II)",
  "grade_groupG3-G4" = "Grau G3-G4 (vs G1-G2)",
  "sexmale" = "Sexo masculino (vs feminino)",
  "age" = "Idade (por ano)",
  "afp_high>=400" = "AFP >=400 (vs <400)",
  "fibrosis_advanced>=5" = "Fibrose Ishak >=5 (vs <5)"
)
forest_df <- cox_multi[cox_multi$term != "(Intercept)", ]
forest_df$rotulo <- ifelse(forest_df$term %in% names(rotulos),
                           rotulos[forest_df$term], forest_df$term)
forest_df$modelo <- factor(forest_df$modelo,
                           levels = c("Modelo A: estágio + grau + sexo + idade",
                                      "Modelo B: A + AFP + fibrose"))
p_forest <- ggplot(forest_df, aes(x = estimate, y = reorder(rotulo, estimate))) +
  geom_point(size = 2.2, colour = "#C0392B") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18,
                 colour = "#C0392B") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ modelo, scales = "free_y") +
  scale_x_log10() +
  labs(title = "Cox multivariado — hazard ratios (IC 95%)",
       subtitle = "Modelo A: n = 342; Modelo B: n = 180 (subgrupo com AFP e fibrose)",
       x = "Hazard ratio (escala log)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(dir_fig, "24_forest_cox_multivariado.png"), p_forest,
       width = 10, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# 6. Correlações clínicas
# -----------------------------------------------------------------------------
vars_cor <- c("age", "afp_log2", "fibrosis_ishak", "stage_num", "grade_num")
labels_cor <- c("Idade", "log2(AFP+1)", "Fibrose Ishak", "Estágio", "Grau")
cor_rho <- matrix(NA_real_, length(vars_cor), length(vars_cor),
                  dimnames = list(labels_cor, labels_cor))
cor_p <- cor_rho
cor_n <- cor_rho
for (i in seq_along(vars_cor)) {
  for (j in seq_along(vars_cor)) {
    x <- tcga[[vars_cor[i]]]; y <- tcga[[vars_cor[j]]]
    ok <- complete.cases(x, y)
    cor_n[i, j] <- sum(ok)
    if (sum(ok) > 5 && sd(x[ok], na.rm = TRUE) > 0 && sd(y[ok], na.rm = TRUE) > 0) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
      cor_rho[i, j] <- ct$estimate
      cor_p[i, j] <- ct$p.value
    }
  }
}
cor_long <- data.frame(
  variavel_1 = rep(labels_cor, each = length(labels_cor)),
  variavel_2 = rep(labels_cor, times = length(labels_cor)),
  rho = as.vector(cor_rho),
  p_valor = as.vector(cor_p),
  n = as.vector(cor_n)
)
write.csv(cor_long, file.path(dir_tab, "19_correlacao_clinica.csv"),
          row.names = FALSE)

png(file.path(dir_fig, "25_correlacao_clinica.png"),
    width = 7, height = 6, units = "in", res = 300)
pheatmap(cor_rho,
         display_numbers = matrix(ifelse(is.na(cor_rho), "",
                                  sprintf("%.2f", cor_rho)), nrow(cor_rho)),
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Correlação de Spearman — fenótipos clínicos (TCGA tumor)",
         fontsize = 10)
dev.off()

# -----------------------------------------------------------------------------
# 7. Resumo
# -----------------------------------------------------------------------------
resumo <- c(
  "FENÓTIPOS AVANÇADOS — RESUMO",
  "============================",
  paste("Tumores primários TCGA:", nrow(tcga)),
  paste("Sobrevida global válida (OS_time > 0):", nrow(surv_df),
        "| eventos:", sum(surv_df$OS)),
  paste("Modelo A (estágio+grau+sexo+idade): n =", nrow(dA),
        "| eventos =", sum(dA$OS)),
  paste("Modelo B (+AFP+fibrose): n =", nrow(dB),
        "| eventos =", sum(dB$OS)),
  "",
  "Missingness (top 8):",
  paste(capture.output(print(head(miss, 8), row.names = FALSE)), collapse = "\n"),
  "",
  "Testes com p < 0,05 (top 15):",
  paste(capture.output(print(head(testes_df[testes_df$p_valor < 0.05, ], 15),
                             row.names = FALSE)), collapse = "\n"),
  "",
  "Log-rank:",
  paste(capture.output(print(logrank_df, row.names = FALSE)), collapse = "\n"),
  "",
  "Cox multivariado:",
  paste(capture.output(print(cox_multi[, c("modelo", "term", "estimate",
                                           "conf.low", "conf.high", "p.value")],
                             row.names = FALSE)), collapse = "\n"),
  "",
  "Sessão R:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(resumo, file.path(dir_qc, "resumo_fenotipos_avancados.txt"))
cat(paste(resumo, collapse = "\n"), "\n")
