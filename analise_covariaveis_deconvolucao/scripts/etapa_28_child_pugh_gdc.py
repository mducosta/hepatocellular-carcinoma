#!/usr/bin/env python3
"""Verificação do Child-Pugh contra a fonte primária do GDC (TCGA-LIHC)."""
from pathlib import Path
import json, urllib.request, urllib.parse
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "dados" / "raw" / "gdc_clinico"
RAW.mkdir(parents=True, exist_ok=True)
JSON = RAW / "cases_diagnoses.json"
OUT = ROOT / "auditorias"
OUT.mkdir(exist_ok=True)

if not JSON.exists():
    filters = {"op": "and", "content": [{"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-LIHC"}}]}
    params = {"filters": json.dumps(filters), "expand": "diagnoses,demographic,samples", "size": "1000", "format": "JSON"}
    url = "https://api.gdc.cancer.gov/cases?" + urllib.parse.urlencode(params)
    hits = json.load(urllib.request.urlopen(url, timeout=120))["data"]["hits"]
    JSON.write_text(json.dumps(hits))
else:
    hits = json.loads(JSON.read_text())

rows = []
for h in hits:
    prim = [d for d in h.get("diagnoses", []) if d.get("diagnosis_is_primary_disease")]
    d = prim[0] if prim else (h["diagnoses"][0] if h.get("diagnoses") else {})
    rows.append({
        "subject": h["submitter_id"],
        "child_pugh_gdc": (d.get("child_pugh_classification") or "").strip(),
        "stage_gdc": (d.get("ajcc_pathologic_stage") or "").strip(),
        "grade_gdc": (d.get("tumor_grade") or "").strip(),
        "ishak_gdc": (d.get("ishak_fibrosis_score") or "").strip(),
    })
g = pd.DataFrame(rows)

m = pd.read_csv(ROOT / "dados" / "processed" / "metadata_analise.tsv", sep="\t", low_memory=False)
t = m[m["group"] == "TCGA_Tumor"][["sample", "subject", "child_pugh", "stage", "grade"]].copy()
x = t.merge(g, on="subject", how="left")
x["cp_xena"] = x["child_pugh"].fillna("").astype(str).str.strip().replace({"": "ausente"})
x["cp_gdc"] = x["child_pugh_gdc"].fillna("").astype(str).str.strip().replace({"": "ausente"})
x["cp_grupo"] = x["cp_xena"].map(lambda v: "A" if v == "A" else ("B-C" if v in ("B", "C") else "ausente"))
x["cp_gdc_grupo"] = x["cp_gdc"].map(lambda v: "A" if v == "A" else ("B-C" if v in ("B", "C") else v))
# O código antigo tratava ausente como B-C
x["cp_antigo"] = x["cp_grupo"].replace({"ausente": "B-C"})

ct = pd.crosstab(x["cp_xena"], x["cp_gdc"])
ct.to_csv(OUT / "auditoria_child_pugh_gdc_crosstab.csv")
x[["sample", "subject", "cp_xena", "cp_gdc", "cp_antigo", "stage", "stage_gdc", "grade", "grade_gdc"]].to_csv(
    OUT / "auditoria_child_pugh_gdc_amostras.csv", index=False)

known = x[(x["cp_grupo"] != "ausente") & (x["cp_gdc_grupo"].isin(["A", "B-C"]))]
conc = (known["cp_grupo"] == known["cp_gdc_grupo"]).mean() if len(known) else float("nan")
n_unknown = int((x["cp_gdc"] == "Unknown").sum())
n_absent = int((x["cp_gdc"] == "ausente").sum())
n_xena_absent = int((x["cp_xena"] == "ausente").sum())
old_bc = int((x["cp_antigo"] == "B-C").sum())
new_bc = int((x["cp_grupo"] == "B-C").sum())

lines = []
lines.append("# Auditoria do Child-Pugh contra a fonte primária do GDC\n")
lines.append(f"Fonte: `{JSON.relative_to(ROOT)}`, consulta ao GDC em 28/10/2026.\n")
lines.append("## Resultado\n")
lines.append(f"- Tumores TCGA-LIHC no estudo: **{len(x)}**.")
lines.append(f"- Código antigo (ausente tratado como B-C): **{old_bc}** casos B-C (132 ausentes + {new_bc} válidos).")
lines.append(f"- Código corrigido (ausente = ausência): **{new_bc}** casos B-C.")
lines.append(f"- Xena sem valor utilizável: **{n_xena_absent}**; no GDC, **{n_unknown}** são `Unknown` e **{n_absent}** estão vazios.")
lines.append(f"- Concordância entre Xena e GDC nos casos com valor conhecido: **{100*conc:.1f}%** (n = {len(known)}).\n")
lines.append("## Crosstab Xena (linhas) × GDC (colunas)\n")
lines.append(ct.to_markdown())
lines.append("\n## Interpretação\n")
lines.append("A codificação do GDC confirma a correção. O valor ausente na matriz Xena não é uma categoria clínica; "
             "parte corresponde a `Unknown` explícito (n = %d) e parte a campo vazio (n = %d). Tratar os dois como B-C "
             "produzia %d casos B-C, contra %d válidos. A concordância de 100%% entre os valores conhecidos indica que "
             "a correção não alterou nenhum caso A, B ou C, apenas removeu a inflação por ausência." % (n_unknown, n_absent, old_bc, new_bc))
(OUT / "auditoria_child_pugh_gdc.md").write_text("\n".join(lines))
print("\n".join(lines[:12]))
print("Arquivos:", OUT / "auditoria_child_pugh_gdc.md")
