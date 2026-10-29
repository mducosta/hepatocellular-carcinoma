#!/usr/bin/env python3
"""Monta a matriz de contagens STAR (GDC) para as 419 amostras do estudo."""
from pathlib import Path
import pandas as pd, glob, os
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "dados" / "raw" / "gdc_counts"
OUT = ROOT / "dados" / "processed" / "gdc_counts"
OUT.mkdir(parents=True, exist_ok=True)
files = sorted(RAW.glob("TCGA-*.tsv"))
print("arquivos:", len(files))
parts = []
for f in files:
    d = pd.read_csv(f, sep="\t", comment="#", usecols=["gene_id", "gene_name", "gene_type", "unstranded"])
    d = d[d["gene_type"] == "protein_coding"]
    d["sample"] = f.stem
    parts.append(d[["gene_id", "gene_name", "sample", "unstranded"]])
m = pd.concat(parts, ignore_index=True)
print("linhas longas:", m.shape)
anot = m[["gene_id", "gene_name"]].drop_duplicates("gene_id").set_index("gene_id")
mat = m.pivot_table(index="gene_id", columns="sample", values="unstranded", aggfunc="sum", fill_value=0)
mat = mat.astype(int)
mat = mat[mat.sum(axis=1) >= 10]
anot = anot.loc[mat.index]
print("matriz:", mat.shape, "genes expressos")
mat.to_csv(OUT / "gdc_counts_proteincoding.tsv.gz", sep="\t", compression="gzip")
anot.to_csv(OUT / "gdc_counts_gene_annotation.tsv", sep="\t")
print("salvo em", OUT)
