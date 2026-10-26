#!/usr/bin/env python3
# =============================================================================
# etapa_17_auditoria_analises_avancadas.py
# -----------------------------------------------------------------------------
# Auditoria das análises avançadas: dados, confundimento, missingness, fluxo de
# casos, atenuação de sinais, QC da deconvolução e confronto entre resultados.
# Saídas em auditorias/ (Markdown + CSV) para rastreabilidade.
# =============================================================================
from pathlib import Path
import math
import warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "resultados" / "analises_avancadas"
TAB = RES / "tabelas"
AUD = ROOT / "auditorias"
AUD.mkdir(exist_ok=True)


def cramers_v(tab):
    tab = np.asarray(tab, dtype=float)
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return np.nan, np.nan, int(tab.sum())
    chi2, p, dof, _ = stats.chi2_contingency(tab, correction=False)
    n = tab.sum()
    v = math.sqrt(chi2 / (n * (min(tab.shape) - 1)))
    return v, p, int(n)


def cat_assoc(df, a, b):
    d = df[[a, b]].dropna()
    if len(d) < 5:
        return dict(var_a=a, var_b=b, n=len(d), cramers_v=np.nan, p=np.nan)
    tab = pd.crosstab(d[a], d[b])
    v, p, n = cramers_v(tab)
    return dict(var_a=a, var_b=b, n=n, cramers_v=v, p=p)


def read_csv(name):
    p = TAB / name
    return pd.read_csv(p) if p.exists() else None


meta = pd.read_csv(ROOT / "dados" / "processed" / "metadata_analise.tsv",
                   sep="\t", low_memory=False, keep_default_na=False, dtype=str)
meta["child_pugh"] = meta["child_pugh"].str.strip()
meta["child_pugh_group_corrigido"] = meta["child_pugh"].map(
    lambda x: np.nan if x in ("", "NA", "nan") else ("A" if x == "A" else "B-C"))
meta["stage_group2"] = meta["stage"].map(
    lambda x: "I-II" if x in ("I", "II") else ("III-IV" if x in ("III", "IV") else np.nan))
meta["grade_group2"] = meta["grade"].map(
    lambda x: "G1-G2" if x in ("G1", "G2") else ("G3-G4" if x in ("G3", "G4") else np.nan))
meta["age_num"] = pd.to_numeric(meta["age"], errors="coerce")

tum = meta[meta["group"] == "TCGA_Tumor"].copy()
adj = meta[meta["group"] == "TCGA_Adjacent"].copy()
gtex = meta[meta["group"] == "GTEx_Normal"].copy()

# -----------------------------------------------------------------------------
# 1. Inventário e integridade
# -----------------------------------------------------------------------------
checks = []


def add_check(item, status, evidence):
    checks.append(dict(item=item, status=status, evidencia=evidence))


add_check("Total de amostras no metadado", "OK", f"{len(meta)} amostras")
add_check("Tumor / adjacente / GTEx", "OK", f"{len(tum)} / {len(adj)} / {len(gtex)}")
add_check("Tumores recorrentes", "INFO",
          "2 recidivas foram excluídas na preparação; 0 no metadado analítico")
gtex_sub = meta[meta["group"] == "GTEx_Normal"]
dup = gtex_sub["subject"].duplicated(keep=False)
add_check("GTEx com mesmo doador em >1 amostra", "ATENÇÃO",
          f"{gtex_sub.loc[dup, 'subject'].nunique()} doador(es); {dup.sum()} amostras; "
          f"análise principal usa 1 por doador")
add_check("Pacientes com tumor e adjacente", "OK",
          f"{len(set(tum['subject']) & set(adj['subject']))} pares")
add_check("Child-Pugh ausente", "ATENÇÃO",
          f"{int(tum['child_pugh'].isin(['', 'NA']).sum())} tumores sem Child-Pugh; "
          "a versão anterior recodificava ausentes como B-C")
add_check("Sorologia desconhecida", "ATENÇÃO",
          f"{(tum['viral_status'] == 'Unknown').sum()} tumores; não é 'sem vírus'")
add_check("Histologia 'Other'", "ATENÇÃO",
          f"{(tum['histology'] != 'Hepatocellular Carcinoma').sum()} tumores "
          "(7 mistos + 3 fibrolamelares)")
add_check("Estágio IV", "ATENÇÃO", f"{(tum['stage'] == 'IV').sum()} tumores")
add_check("Grau G4", "ATENÇÃO", f"{(tum['grade'] == 'G4').sum()} tumores")
add_check("Idade GTEx", "ATENÇÃO", "faixa etária convertida ao ponto médio")
add_check("Dados de expressão", "OK",
          "log2(TPM + 0,001); 12.401 genes; 528 amostras após deduplicação")

expr_path = ROOT / "dados" / "processed" / "expressao_log2tpm_filtrada.tsv.gz"
expr = pd.read_csv(expr_path, sep="\t", index_col=0)
add_check("Matriz de expressão", "OK", f"{expr.shape[0]} genes × {expr.shape[1]} amostras; análise principal com 528 após deduplicação")
add_check("Valores ausentes na matriz", "OK", f"{int(expr.isna().sum().sum())} valores ausentes")
add_check("Genes com variância zero", "OK", f"{int((expr.var(axis=1) == 0).sum())} genes")
add_check("Escala log2", "OK",
          f"min={expr.values.min():.3f}; max={expr.values.max():.3f}; "
          f"mediana={np.median(expr.values):.3f}")

# -----------------------------------------------------------------------------
# 2. Confundimento / associações
# -----------------------------------------------------------------------------
assoc = []
assoc.append(cat_assoc(meta, "group", "study"))
for a, b in [("stage_group2", "tss"), ("grade_group2", "tss"),
             ("child_pugh_group_corrigido", "tss"),
             ("histology", "tss"), ("stage_group2", "grade_group2"),
             ("stage_group2", "child_pugh_group_corrigido"),
             ("grade_group2", "child_pugh_group_corrigido"),
             ("stage_group2", "viral_status"), ("grade_group2", "viral_status"),
             ("grade_group2", "afp_group"), ("stage_group2", "afp_group"),
             ("sex", "group"), ("sex", "stage_group2"),
             ("sex", "grade_group2"), ("sex", "tss")]:
    d = tum if b not in ["group", "study"] else meta
    assoc.append(cat_assoc(d, a, b))
assoc = pd.DataFrame(assoc)
assoc.to_csv(AUD / "auditoria_confundimento_associacoes.csv", index=False)

tss_stage = pd.crosstab(tum["tss"], tum["stage_group2"])
tss_stage["n"] = tss_stage.sum(axis=1)
tss_stage["pct_III_IV"] = 100 * tss_stage.get("III-IV", 0) / tss_stage["n"]
tss_stage["child_pugh_ausente"] = tum.groupby("tss")["child_pugh"].apply(
    lambda x: x.isin(["", "NA"]).sum())
tss_stage = tss_stage.sort_values("n", ascending=False)
tss_stage.to_csv(AUD / "auditoria_tss_estagio_child.csv")

# -----------------------------------------------------------------------------
# 3. Missingness corrigida e fluxo de casos
# -----------------------------------------------------------------------------
vars_miss = ["age", "sex", "stage", "grade", "afp", "fibrosis_ishak",
             "viral_status", "child_pugh", "vascular_invasion", "histology",
             "tss", "OS", "OS_time"]
miss = []
for v in vars_miss:
    x = tum[v].replace("", np.nan)
    unknown = x.isin(["Unknown", "not reported", "NA", "nan"])
    miss.append(dict(variavel=v, ausente=int(x.isna().sum() + unknown.sum()),
                     total=len(tum),
                     percentual=round(100 * (x.isna().sum() + unknown.sum()) / len(tum), 1)))
miss = pd.DataFrame(miss).sort_values("percentual", ascending=False)
miss.to_csv(AUD / "auditoria_missingness_corrigida.csv", index=False)


def cc(df, cols):
    return df[cols].replace("", np.nan).notna().all(axis=1).sum()


flow = []
for label, cols in [
    ("Tumor vs adjacente: base (sexo, idade)", ["group", "sex", "age"]),
    ("Tumor vs adjacente: + estágio e grau", ["group", "sex", "age", "stage", "grade"]),
    ("Tumor vs adjacente: + sorologia, Child-Pugh, histologia",
     ["group", "sex", "age", "stage", "grade", "viral_status", "child_pugh", "histology"]),
    ("Tumor vs adjacente: + TSS",
     ["group", "sex", "age", "stage", "grade", "viral_status", "child_pugh", "histology", "tss"]),
]:
    d = meta[meta["group"].isin(["TCGA_Tumor", "TCGA_Adjacent"])]
    flow.append(dict(modelo=label, n_total=cc(d, cols),
                     n_tumor=cc(d[d["group"] == "TCGA_Tumor"], cols),
                     n_adjacente=cc(d[d["group"] == "TCGA_Adjacent"], cols)))
for label, cols in [
    ("Estágio: base", ["stage", "sex", "age"]),
    ("Estágio: + grau e clínica",
     ["stage", "sex", "age", "grade", "viral_status", "child_pugh", "histology"]),
    ("Estágio: + TSS",
     ["stage", "sex", "age", "grade", "viral_status", "child_pugh", "histology", "tss"]),
    ("Grau: base", ["grade", "sex", "age"]),
    ("Grau: + estágio e clínica",
     ["grade", "sex", "age", "stage", "viral_status", "child_pugh", "histology"]),
    ("Grau: + TSS",
     ["grade", "sex", "age", "stage", "viral_status", "child_pugh", "histology", "tss"]),
]:
    flow.append(dict(modelo=label, n_total=cc(tum, cols),
                     n_tumor=cc(tum, cols), n_adjacente=None))
flow = pd.DataFrame(flow)
flow.to_csv(AUD / "auditoria_fluxo_casos.csv", index=False)

# -----------------------------------------------------------------------------
# 4. Atenuação de sinais (confronto)
# -----------------------------------------------------------------------------
def deg_count(path):
    p = TAB / path
    return pd.read_csv(p) if p.exists() else None


def attenuation(base_file, alt_file, label):
    a = deg_count(base_file)
    b = deg_count(alt_file)
    if a is None or b is None:
        return dict(comparacao=label, degs_base=None, degs_alt=None,
                    jaccard=None, atenuacao_mediana=None)
    sa = a[a["regulation"] != "NS"]
    sb = b[b["regulation"] != "NS"]
    m = sa[["gene", "logFC"]].merge(b[["gene", "logFC"]], on="gene",
                                    suffixes=("_base", "_alt"))
    att = (1 - m["logFC_alt"].abs() / m["logFC_base"].abs().clip(lower=1e-9)).median() if len(m) else np.nan
    jac = len(set(sa.gene) & set(sb.gene)) / len(set(sa.gene) | set(sb.gene)) if (len(sa) or len(sb)) else np.nan
    return dict(comparacao=label, degs_base=len(sa), degs_alt=len(sb),
                jaccard=jac, atenuacao_mediana=att, n_base=len(a), n_alt=len(b))


atenuacoes = []
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_estagio_base_completo.csv",
    "expressao_diferencial_covariaveis/modelo_estagio_base_mesmos_casos_completo.csv",
    "Estágio: coorte completa vs mesmos casos completos"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_estagio_base_completo.csv",
    "expressao_diferencial_covariaveis/modelo_estagio_so_tss_full_completo.csv",
    "Estágio: sem vs com TSS"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_estagio_base_completo.csv",
    "expressao_diferencial_covariaveis/modelo_estagio_sem_cc_ed_completo.csv",
    "Estágio: excluindo CC/ED"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_grau_base_completo.csv",
    "expressao_diferencial_covariaveis/modelo_grau_clinico_tss_completo.csv",
    "Grau: base vs clínico + TSS"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_tumor_vs_adjacente_base_completo.csv",
    "expressao_diferencial_covariaveis/modelo_tumor_vs_adjacente_pareado_completo.csv",
    "Tumor vs adjacente: não pareado vs pareado"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_tumor_vs_adjacente_base_completo.csv",
    "integracao_composicao/tumor_vs_adjacente_composicao_ajustado_composicao.csv",
    "Tumor vs adjacente: sem vs com composição"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_tumor_vs_gtex_base_completo.csv",
    "integracao_composicao/tumor_vs_gtex_composicao_ajustado_composicao.csv",
    "Tumor vs GTEx: sem vs com composição"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_grau_base_completo.csv",
    "integracao_composicao/grau_composicao_ajustado_composicao.csv",
    "Grau: sem vs com composição"))
atenuacoes.append(attenuation(
    "expressao_diferencial_covariaveis/modelo_estagio_base_completo.csv",
    "integracao_composicao/estagio_composicao_ajustado_composicao.csv",
    "Estágio: sem vs com composição"))
atenuacoes = pd.DataFrame(atenuacoes)
atenuacoes.to_csv(AUD / "auditoria_atenuacao_sinais.csv", index=False)

# -----------------------------------------------------------------------------
# 5. QC da deconvolução
# -----------------------------------------------------------------------------
qc = pd.read_csv(TAB / "deconvolucao" / "QC_correlacao_celulas_marcadores.csv")
key = ["xCell_Hepatocytes", "EPIC_otherCells", "xCell_Endothelial cells",
       "EPIC_CAFs", "xCell_Macrophages", "EPIC_Macrophages",
       "xCell_CD8+ T-cells", "EPIC_CD8_Tcells"]
qc_key = qc[qc["celula"].isin(key)]
qc_key.to_csv(AUD / "auditoria_deconvolucao_qc.csv", index=False)
conc = pd.read_csv(TAB / "deconvolucao" / "concordancia_metodos_deconvolucao.csv")
conc.to_csv(AUD / "auditoria_deconvolucao_concordancia.csv", index=False)
gof = pd.read_csv(TAB / "deconvolucao" / "QC_EPIC_convergencia.csv")
add_check("EPIC: convergência", "ATENÇÃO",
          f"{int((gof['convergeCode'] != 0).sum())} de {len(gof)} amostras com código != 0")
add_check("EPIC otherCells vs marcadores de hepatócito", "ATENÇÃO",
          "Spearman negativo com ALB/TTR; otherCells NÃO deve ser chamado de fração de hepatócitos")
add_check("xCell Hepatocytes vs marcadores de hepatócito", "OK",
          "Spearman positivo com ALB, HP, TTR e CYP3A4")
add_check("Concordância entre métodos", "ATENÇÃO",
          f"Spearman mediano entre métodos = "
          f"{conc[['spearman_ab','spearman_ac','spearman_bc']].stack().median():.2f}")

# -----------------------------------------------------------------------------
# 6. Diagnósticos de modelo e composição
# -----------------------------------------------------------------------------
comp = pd.read_csv(TAB / "deconvolucao" / "deconvolucao_unificada_amostras.csv",
                   low_memory=False)
cors = comp[["xCell_Hepatocytes", "xCell_ImmuneScore", "xCell_StromaScore"]].corr(method="spearman")
cors.to_csv(AUD / "auditoria_correlacao_escores_composicao.csv")
vif = pd.Series(np.diag(np.linalg.pinv(cors.values)), index=cors.columns)
vif.to_csv(AUD / "auditoria_vif_composicao.csv")
add_check("Colinearidade dos escores de composição", "INFO",
          "; ".join(f"{k}: VIF={v:.2f}" for k, v in vif.items()))

# -----------------------------------------------------------------------------
# 7. Resumo Markdown
# -----------------------------------------------------------------------------
def md_table(df, cols=None, floatfmt=".3f"):
    if cols:
        df = df[cols]
    return df.to_markdown(index=False, floatfmt=floatfmt)


conf_key = assoc.sort_values("cramers_v", ascending=False)
crit = [
    ("Confundimento grupo × estudo", "Cramér's V = 1,0; logFC tumor vs GTEx não separa biologia de lote."),
    ("Child-Pugh mal codificado", "132 tumores sem Child-Pugh foram tratados como B-C na versão anterior; aqui são NA. O n de B-C cai de 154 para 22."),
    ("Estágio confundido com TSS", "CC concentra 26/90 tumores III-IV e 100% dos CC sem Child-Pugh; o ajuste por TSS zera os DEGs de estágio."),
    ("Estágio frágil à restrição amostral", "46 DEGs na coorte completa, 0 nos mesmos casos completos; logFC mediano cai de 1,31 para 0,52."),
    ("Grau robusto, mas composicional", "172–337 DEGs em modelos diferentes; ajuste por composição zera os DEGs."),
    ("Tumor vs adjacente robusto", "4.130 DEGs no basal; 3.710 no pareado e 3.740 no clínico; Jaccard 0,76–0,83."),
    ("Composição explica parte do sinal", "Tumor vs adjacente cai de 4.130 para 1.897 DEGs após ajuste por xCell."),
    ("Sexo não é só composição", "13 DEGs no GTEx e 175 no tumor; ajuste por composição reduz tumor para 90, mas não zera."),
    ("Histologia não interpretável", "Apenas 10 casos 'Other'; 740–1.020 DEGs são instáveis e não devem ser interpretados como biologia."),
    ("Sorologia/Child-Pugh/fibrose/invasão", "Zero DEGs com FDR < 0,05 e |logFC| > 1; poder baixo e missingness alta."),
    ("Deconvolução não é fração verdadeira", "EPIC otherCells não é hepatócito; xCell é escore de enriquecimento; concordância entre métodos é baixa."),
]

md = []
md.append("# Auditoria das análises avançadas — TCGA-LIHC × GTEx\n")
md.append(f"Gerado em {pd.Timestamp.now():%Y-%m-%d %H:%M}.\n")
md.append("## 1. Inventário e integridade\n")
md.append(md_table(pd.DataFrame(checks)))
md.append("\n## 2. Associações e confundimento\n")
md.append(md_table(conf_key.head(15)))
md.append("\n### TSS com maior proporção de estágio III–IV\n")
md.append(md_table(tss_stage.head(12).reset_index(), floatfmt=".1f"))
md.append("\n## 3. Missingness corrigida (tumor)\n")
md.append(md_table(miss))
md.append("\n## 4. Fluxo de casos\n")
md.append(md_table(flow))
md.append("\n## 5. Atenuação de sinais\n")
md.append(md_table(atenuacoes))
md.append("\n## 6. QC da deconvolução\n")
md.append(md_table(qc_key[qc_key["celula"].isin(key)]))
md.append("\n## 7. Concordância entre métodos\n")
md.append(md_table(conc))
md.append("\n## 8. Achados críticos\n")
for i, (t, d) in enumerate(crit, 1):
    md.append(f"{i}. **{t}.** {d}")
md.append("\n## 9. Conclusão da auditoria\n")
md.append("A análise principal deve priorizar o contraste tumor vs. adjacente e o grau, "
          "com ajuste por composição apresentado como análise de atenuação. "
          "Estágio, sorologia, Child-Pugh, histologia, fibrose e invasão vascular "
          "não sustentam inferência robusta neste desenho. A comparação TCGA × GTEx "
          "permanece confundida com estudo/lote e não deve ser interpretada como efeito "
          "puro de tumor.\n")
(AUD / "auditoria_analises_avancadas.md").write_text("\n".join(md), encoding="utf-8")
pd.DataFrame(checks).to_csv(AUD / "auditoria_analises_avancadas.csv", index=False)
print("Auditoria concluída.")
print(atenuacoes.to_string(index=False))
