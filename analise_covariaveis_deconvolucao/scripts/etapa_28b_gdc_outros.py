#!/usr/bin/env python3
"""Verifica estágio, grau e fibrose contra a fonte primária do GDC."""
from pathlib import Path
import json, re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "auditorias"
OUT.mkdir(exist_ok=True)
hits = json.loads((ROOT / "dados" / "raw" / "gdc_clinico" / "cases_diagnoses.json").read_text())

rows = []
for h in hits:
    prim = [d for d in h.get("diagnoses", []) if d.get("diagnosis_is_primary_disease")]
    d = prim[0] if prim else (h["diagnoses"][0] if h.get("diagnoses") else {})
    stage = (d.get("ajcc_pathologic_stage") or "").strip()
    stage_g = None
    if re.search(r"Stage (I|II)$", stage): stage_g = "I-II"
    elif re.search(r"Stage (III|IV)", stage): stage_g = "III-IV"
    grade = (d.get("tumor_grade") or "").strip()
    ishak = (d.get("ishak_fibrosis_score") or "").strip()
    ishak_n = int(re.match(r"(\d+)", ishak).group(1)) if re.match(r"(\d+)", ishak) else None
    rows.append(dict(subject=h["submitter_id"], stage_gdc=stage_g, stage_raw=stage,
                     grade_gdc=grade, ishak_gdc=ishak_n))
g = pd.DataFrame(rows)

m = pd.read_csv(ROOT / "dados" / "processed" / "metadata_analise.tsv", sep="\t", low_memory=False)
t = m[m["group"] == "TCGA_Tumor"].copy()
t["stage_g"] = t["stage"].map(lambda s: "I-II" if s in ("I", "II") else ("III-IV" if s in ("III", "IV") else None))
t["grade_g"] = t["grade"].map(lambda s: s if s in ("G1", "G2", "G3", "G4") else None)
t["ishak_adv"] = t["fibrosis_ishak"].map(lambda x: None if pd.isna(x) else (">=5" if x >= 5 else "<5"))
g["ishak_gdc_adv"] = g["ishak_gdc"].map(lambda x: None if pd.isna(x) else (">=5" if x >= 5 else "<5"))

x = t[["sample", "subject", "stage_g", "grade_g", "ishak_adv"]].merge(g, on="subject", how="left")
res = []
for nome, a, b in [("Estágio", "stage_g", "stage_gdc"), ("Grau", "grade_g", "grade_gdc"),
                   ("Fibrose", "ishak_adv", "ishak_gdc_adv")]:
    d = x[[a, b]].dropna()
    if len(d):
        conc = (d[a] == d[b]).mean()
        res.append(dict(variavel=nome, n_comparado=len(d), concordancia=round(100*conc, 1),
                        xena_ausente=int(x[a].isna().sum()), gdc_ausente=int(x[b].isna().sum())))
        pd.crosstab(d[a], d[b]).to_csv(OUT / f"auditoria_gdc_{nome.lower()}_crosstab.csv")
r = pd.DataFrame(res)
r.to_csv(OUT / "auditoria_gdc_outros.csv", index=False)
lines = ["# Auditoria de estágio, grau e fibrose contra o GDC\n",
         "Fonte: `dados/raw/gdc_clinico/cases_diagnoses.json`, diagnóstico primário.\n",
         "| Variável | n comparado | Concordância (%) | Xena ausente | GDC ausente |",
         "|---|---:|---:|---:|---:|"]
for _, z in r.iterrows():
    lines.append(f"| {z['variavel']} | {int(z['n_comparado'])} | {z['concordancia']} | {int(z['xena_ausente'])} | {int(z['gdc_ausente'])} |")
(OUT / "auditoria_gdc_outros.md").write_text("\n".join(lines))
print(r.to_string(index=False))
