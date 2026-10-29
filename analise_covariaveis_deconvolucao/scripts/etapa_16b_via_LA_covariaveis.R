# =============================================================================
# etapa_16b_via_LA_covariaveis.R
# -----------------------------------------------------------------------------
# Reanálise da via KEGG hsa05417 (Lipídios e Aterosclerose) com as mesmas
# covariáveis e a mesma deconvolução usadas no artigo, para conectar o achado
# original do seminário à análise de robustez.
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
SCRIPT_DIR <- if (length(script_arg) == 1) dirname(normalizePath(script_arg)) else getwd()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "utilitarios_analises_avancadas.R"))
dir_proc <- file.path(ROOT, "dados", "processed")
dir_out <- file.path(ROOT, "resultados", "analises_avancadas")
dir_tab <- file.path(dir_out, "tabelas", "via_LA")
dir_qc <- file.path(dir_out, "qc")
dir.create(dir_tab, recursive=TRUE, showWarnings=FALSE)
X <- carregar_expressao(file.path(dir_proc, "expressao_log2tpm_filtrada.tsv.gz"))
la <- fread(file.path(dir_proc, "expressao_LA_pathway.tsv.gz"), sep="\t", data.table=FALSE)
la_genes <- intersect(la$gene, rownames(X))
X <- X[la_genes, , drop=FALSE]
uni <- fread(file.path(dir_out, "tabelas", "deconvolucao", "deconvolucao_unificada_amostras.csv"), sep=",", data.table=FALSE)
uni <- uni[uni$sample %in% colnames(X),]
for (v in c("group","sex","stage_group","grade_group","viral_status","child_pugh_group","histology_group","tss_collapsed")) uni[[v]] <- factor(uni[[v]])
for (v in c("xCell_Hepatocytes","xCell_ImmuneScore","xCell_StromaScore")) uni[[paste0("z_",v)]] <- as.numeric(scale(uni[[v]]))
comp <- c("z_xCell_Hepatocytes","z_xCell_ImmuneScore","z_xCell_StromaScore")
d_tumor <- uni[uni$group=="TCGA_Tumor",]; d_adj <- uni[uni$group=="TCGA_Adjacent",]; d_gtex <- uni[uni$group=="GTEx_Normal",]
d_tcga <- uni[uni$group %in% c("TCGA_Tumor","TCGA_Adjacent"),]; d_tg <- uni[uni$group %in% c("TCGA_Tumor","GTEx_Normal"),]
models <- list(
 list(nome="LA_tumor_vs_adjacente_base", data=d_tcga, gc="group", lev=c("TCGA_Adjacent","TCGA_Tumor"), cov=c("sex","age")),
 list(nome="LA_tumor_vs_adjacente_clinico", data=d_tcga, gc="group", lev=c("TCGA_Adjacent","TCGA_Tumor"), cov=c("sex","age","stage_group","grade_group","viral_status","child_pugh_group","histology_group")),
 list(nome="LA_tumor_vs_adjacente_composicao", data=d_tcga, gc="group", lev=c("TCGA_Adjacent","TCGA_Tumor"), cov=c("sex","age",comp)),
 list(nome="LA_tumor_vs_gtex_base", data=d_tg, gc="group", lev=c("GTEx_Normal","TCGA_Tumor"), cov=c("sex","age")),
 list(nome="LA_tumor_vs_gtex_composicao", data=d_tg, gc="group", lev=c("GTEx_Normal","TCGA_Tumor"), cov=c("sex","age",comp)),
 list(nome="LA_estagio_clinico", data=d_tumor, gc="stage_group", lev=c("I-II","III-IV"), cov=c("sex","age","grade_group","viral_status","child_pugh_group","histology_group")),
 list(nome="LA_grau_clinico", data=d_tumor, gc="grade_group", lev=c("G1-G2","G3-G4"), cov=c("sex","age","stage_group","viral_status","child_pugh_group","histology_group")),
 list(nome="LA_sexo_gtex", data=d_gtex, gc="sex", lev=c("female","male"), cov=c("age")),
 list(nome="LA_sexo_tumor", data=d_tumor, gc="sex", lev=c("female","male"), cov=c("age","stage_group","grade_group","viral_status","histology_group","tss_collapsed"))
)
res <- list(); tops <- list()
for (m in models) {
  out <- do.call(rodar_limma, c(list(X=X, d=m$data, group_col=m$gc, levels=m$lev, covariates=m$cov), list(fdr_cut=0.05,lfc_cut=1)))
  if (is.null(out)) next
  out$modelo <- m$nome
  fwrite(out, file.path(dir_tab, paste0(m$nome,"_completo.csv")))
  sig <- out[out$regulation!="NS",]
  fwrite(sig, file.path(dir_tab, paste0(m$nome,"_significativos.csv")))
  res[[m$nome]] <- data.frame(modelo=m$nome, n=length(unique(out$gene)), up=sum(out$regulation=="Up"), down=sum(out$regulation=="Down"), total=sum(out$regulation!="NS"), min_fdr=min(out$adj.P.Val))
  tops[[m$nome]] <- head(out[out$regulation!="NS",][order(out$adj.P.Val[out$regulation!="NS"]),], 10)
}
res <- do.call(rbind,res); rownames(res)<-NULL
fwrite(res, file.path(dir_tab,"resumo_via_LA_covariaveis.csv"))
writeLines(c("VIA LA — REANÁLISE COM COVARIÁVEIS", paste("Genes:",length(la_genes)), paste0(res$modelo,": ",res$total," DEGs (",res$up," up / ",res$down," down)")), file.path(dir_qc,"log_via_LA_covariaveis.txt"))
print(res)
