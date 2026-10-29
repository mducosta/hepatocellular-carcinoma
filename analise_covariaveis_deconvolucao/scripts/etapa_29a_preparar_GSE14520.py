#!/usr/bin/env python3
"""Prepara GSE14520 (GPL3921) para validação externa."""
from pathlib import Path
import gzip, io, csv, re
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "dados" / "raw" / "externo"
OUT = ROOT / "resultados" / "analises_avancadas" / "tabelas" / "validacao_externa"
OUT.mkdir(parents=True, exist_ok=True)

# --- matriz ---
with gzip.open(RAW / "GSE14520-GPL3921_series_matrix.txt.gz", "rt", encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()
meta_lines = [l for l in lines if l.startswith("!Sample_")]
idx = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
tab = pd.read_csv(io.StringIO("".join(lines[idx+1:end])), sep="\t", quotechar='"')
tab = tab.rename(columns={"ID_REF": "probe"})
print("matriz:", tab.shape)

def campos(prefix):
    for l in meta_lines:
        if l.startswith(prefix):
            return next(csv.reader([l.rstrip("\n")], delimiter="\t"))[1:]
    return []

gsm = campos("!Sample_geo_accession")
title = campos("!Sample_title")
char_lines = [l for l in meta_lines if l.startswith("!Sample_characteristics_ch1")]
chars = [next(csv.reader([l.rstrip()], delimiter=chr(9)))[1:] for l in char_lines[:3]]
# As três linhas de características são Tissue, Disease state e Individual
meta = pd.DataFrame({"gsm": gsm, "title": title,
                     "tissue_raw": chars[0], "disease": chars[1], "individual_raw": chars[2]})
meta["tissue"] = np.select(
    [meta["tissue_raw"].str.contains("Non-Tumor", case=False),
     meta["tissue_raw"].str.contains("Tumor", case=False)],
    ["Non-Tumor", "Tumor"], default="Other")
meta["individual"] = meta["individual_raw"].str.replace("Individual: ", "", regex=False).str.strip()
meta["pair"] = meta["individual"].str[:-1]
print(meta["tissue"].value_counts())

# --- anotação GPL3921 ---
ann = pd.read_csv(RAW / "GPL3921.annot.gz", sep="\t", skiprows=27, dtype=str,
                  usecols=["ID", "Gene symbol"], quotechar='"')
ann = ann.rename(columns={"ID": "probe", "Gene symbol": "gene"})
ann = ann.dropna(subset=["gene"])
ann = ann[ann["gene"].str.strip() != ""]
ann = ann.assign(gene=ann["gene"].str.split("///")).explode("gene")
ann["gene"] = ann["gene"].str.strip()
print("anotação:", ann.shape, "genes únicos:", ann["gene"].nunique())

# --- clínico suplementar ---
sup = pd.read_csv(RAW / "GSE14520_Extra_Supplement.txt.gz", sep="\t")
sup = sup.rename(columns={"Affy_GSM": "gsm", "Gender": "sex", "Age": "age", "TNM staging": "tnm"})
sup = sup[["gsm", "Tissue Type", "sex", "age", "tnm", "HBV viral status", "AFP (>/<=300ng/ml)"]]
meta = meta.merge(sup, on="gsm", how="left")
meta["stage"] = meta["tnm"].astype(str).str.replace(r"^III[A-C]$", "III", regex=True)
meta["stage_group"] = meta["stage"].map(lambda s: "I-II" if s in ("I", "II") else ("III-IV" if s in ("III", "IV") else np.nan))
print("estágio:", meta["stage_group"].value_counts(dropna=False).to_dict())

# --- colapsa sondas por gene (maior expressão média) ---
X = tab.set_index("probe").drop(columns=[c for c in ["IDENTIFIER"] if c in tab.columns], errors="ignore")
X = X.apply(pd.to_numeric, errors="coerce")
X = X.loc[:, [c for c in X.columns if c in meta["gsm"].values]]
ann2 = ann[ann["probe"].isin(X.index)]
mean_expr = X.mean(axis=1)
ann2 = ann2.assign(mean_expr=ann2["probe"].map(mean_expr))
ann2 = ann2.sort_values("mean_expr", ascending=False).drop_duplicates("gene")
Xg = X.loc[ann2["probe"].values, :]
Xg.index = ann2["gene"].values
Xg = Xg.groupby(level=0).max()
print("matriz gênica:", Xg.shape)

meta.to_csv(OUT / "GSE14520_metadata.tsv", sep="\t", index=False)
Xg.to_csv(OUT / "GSE14520_expression.tsv", sep="\t")
print("salvo em", OUT)
