#!/usr/bin/env python3
"""Prepara GSE62232, GSE36376 e GSE76427 para validação em coortes independentes."""
from pathlib import Path
import gzip, io, csv
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "dados" / "raw" / "externo"
OUT = ROOT / "resultados" / "analises_avancadas" / "tabelas" / "validacao_externa"
OUT.mkdir(parents=True, exist_ok=True)

ANNOT = {
    "GPL570": RAW / "GPL570_annotation.tsv",
    "GPL10558": RAW / "GPL10558.annot.gz",
    "GPL3921": RAW / "GPL3921.annot.gz",
}


def campos(lines, prefix):
    for l in lines:
        if l.startswith(prefix):
            return next(csv.reader([l.rstrip("\n")], delimiter="\t"))[1:]
    return []


def ler_series(gse):
    with gzip.open(RAW / f"{gse}_series_matrix.txt.gz", "rt", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    idx = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    tab = pd.read_csv(io.StringIO("".join(lines[idx+1:end])), sep="\t", quotechar='"')
    tab = tab.rename(columns={tab.columns[0]: "probe"})
    meta = pd.DataFrame({"gsm": campos(lines, "!Sample_geo_accession"),
                         "title": campos(lines, "!Sample_title")})
    # todas as linhas de características, na ordem
    chars = []
    for l in lines:
        if l.startswith("!Sample_characteristics_ch1"):
            chars.append(next(csv.reader([l.rstrip("\n")], delimiter="\t"))[1:])
    for i, ch in enumerate(chars):
        meta[f"ch{i}"] = ch
    plat = campos(lines, "!Sample_platform_id")
    meta["platform"] = plat[0] if plat else ""
    return tab, meta


def carregar_anotacao(platform):
    if platform == "GPL570":
        ann = pd.read_csv(ANNOT[platform], sep="\t", dtype=str)
        ann = ann.rename(columns={"ID": "probe", "Gene Symbol": "gene"})
    else:
        f = ANNOT[platform]
        ann = pd.read_csv(f, sep="\t", skiprows=28, dtype=str, usecols=["ID", "Gene symbol"],
                          quotechar='"').rename(columns={"ID": "probe", "Gene symbol": "gene"})
    ann = ann.dropna(subset=["gene"])
    ann = ann[ann["gene"].str.strip() != ""]
    ann = ann.assign(gene=ann["gene"].str.split("///")).explode("gene")
    ann["gene"] = ann["gene"].str.strip()
    return ann


def colapsar(tab, platform, samples):
    ann = carregar_anotacao(platform)
    X = tab.set_index("probe")
    X = X.loc[:, [c for c in X.columns if c in samples]]
    X = X.apply(pd.to_numeric, errors="coerce")
    mean_expr = X.mean(axis=1)
    a = ann[ann["probe"].isin(X.index)].assign(mean_expr=lambda d: d["probe"].map(mean_expr))
    a = a.sort_values("mean_expr", ascending=False).drop_duplicates("gene")
    Xg = X.loc[a["probe"].values, :]
    Xg.index = a["gene"].values
    Xg = Xg.groupby(level=0).max()
    # Plataformas entregues em intensidade linear (GSE62232/GSE76427) precisam de log2
    if np.nanmedian(Xg.values) > 50 or np.nanmax(Xg.values) > 100:
        Xg = np.log2(np.clip(Xg, 0, None) + 1)
    return Xg


def prep(gse, tissue_map, grade_fn, stage_fn, sex_fn, age_fn, pair_fn=None):
    tab, meta = ler_series(gse)
    char_cols = [c for c in meta.columns if c.startswith("ch")]
    recs = []
    for _, r in meta.iterrows():
        d = {}
        for c in char_cols:
            v = str(r[c]).strip().strip(chr(34)) if pd.notna(r[c]) else ""
            if ":" in v:
                k, val = v.split(":", 1)
                d[k.strip().strip(chr(34)).lower()] = val.strip().strip(chr(34))
        recs.append(d)
    meta["_chars"] = recs
    meta["tissue"] = meta["_chars"].map(lambda d: tissue_map.get(str((d or {}).get("tissue", "")).strip().lower()))
    meta["grade"] = meta["_chars"].map(grade_fn)
    meta["stage_raw"] = meta["_chars"].map(stage_fn)
    meta["sex"] = meta["_chars"].map(sex_fn)
    meta["age"] = meta["_chars"].map(age_fn)
    meta["pair"] = meta["_chars"].map(pair_fn) if pair_fn else meta["title"]
    meta = meta[meta["tissue"].isin(["Tumor", "Non-Tumor"])].copy()
    meta = meta[~meta["title"].str.contains(r"\.rep", regex=True)]
    platform = meta["platform"].iloc[0]
    X = colapsar(tab, platform, set(meta["gsm"]))
    meta = meta[meta["gsm"].isin(X.columns)].copy()
    X = X.loc[:, meta["gsm"]]
    meta["stage_group"] = meta["stage_raw"].map(lambda s: "I-II" if str(s).upper() in ("1", "2", "I", "II") else
                                                   ("III-IV" if str(s).upper() in ("3", "3A", "3B", "3C", "4", "III", "IIIA", "IIIB", "IIIC", "IV") else None))
    meta["grade_group"] = meta["grade"].map(lambda s: "G1-G2" if s and str(s).lower() in ("i-ii", "i", "ii", "g1", "g2", "1", "2") else
                                                       ("G3-G4" if s and str(s).lower() in ("iii-iv", "iii", "iv", "g3", "g4", "3", "4") else None))
    out = meta[["gsm", "title", "tissue", "grade", "grade_group", "stage_raw", "stage_group", "sex", "age", "pair"]].copy()
    out.insert(0, "cohort", gse)
    out.to_csv(OUT / f"{gse}_metadata.tsv", sep="	", index=False)
    X.to_csv(OUT / f"{gse}_expression.tsv", sep="	")
    print(gse, X.shape, out["tissue"].value_counts().to_dict(), out["grade_group"].value_counts(dropna=False).to_dict())
    print(out["stage_group"].value_counts(dropna=False).to_dict(), out["sex"].value_counts(dropna=False).to_dict())


prep("GSE62232",
     {"non tumor liver": "Non-Tumor", "hepatocellular carcinoma": "Tumor"},
     lambda d: d.get("edmonson"), lambda d: d.get("ajcc stage") or d.get("tnm staging"),
     lambda d: d.get("gender"), lambda d: d.get("age"))
prep("GSE36376",
     {"adjacent non-tumor liver": "Non-Tumor", "liver tumor": "Tumor"},
     lambda d: d.get("grade") or d.get("edmondson"), lambda d: d.get("ajcc stage"),
     lambda d: d.get("gender") or d.get("sex"), lambda d: d.get("age") or d.get("age (years)"))
prep("GSE76427",
     {"primary hepatocellular carcinoma tumor": "Tumor", "adjacent non-tumor liver tissue": "Non-Tumor"},
     lambda d: d.get("grade") or d.get("edmondson"), lambda d: d.get("tnm_staging_clinical"),
     lambda d: {"1": "M", "2": "F"}.get(d.get("gender (1=m, 2=f)"), d.get("gender")),
     lambda d: d.get("age (years)"), pair_fn=lambda d: d.get("patient id"))
