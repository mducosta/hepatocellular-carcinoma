#!/usr/bin/env python3
"""Baixa as contagens STAR (GDC) para as amostras TCGA-LIHC do estudo."""
import json, os, sys, time, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "dados", "raw", "gdc_counts")
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(os.path.join(ROOT, "dados", "processed", "metadata_analise.tsv"),
                   sep="\t", low_memory=False)
cohort = meta[meta["group"].isin(["TCGA_Tumor", "TCGA_Adjacent"])][["sample", "group", "subject"]].copy()
cohort["key"] = cohort["sample"].str[:15]

filters = {"op": "and", "content": [
    {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-LIHC"}},
    {"op": "=", "content": {"field": "data_category", "value": "Transcriptome Profiling"}},
    {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
    {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
]}
params = {"filters": json.dumps(filters), "fields": "file_id,file_name,file_size,cases.samples.sample_type,cases.samples.submitter_id",
          "size": "2000", "format": "JSON"}
url = "https://api.gdc.cancer.gov/files?" + urllib.parse.urlencode(params)
hits = json.load(urllib.request.urlopen(url, timeout=120))["data"]["hits"]
rows = []
for h in hits:
    st = h["cases"][0]["samples"][0]
    key = st["submitter_id"][:15]
    rows.append({"key": key, "sample_gdc": st["submitter_id"], "sample_type": st["sample_type"],
                 "file_id": h["file_id"], "file_name": h["file_name"], "file_size": h["file_size"]})
man = pd.DataFrame(rows)
man = man[man["sample_type"].isin(["Primary Tumor", "Solid Tissue Normal"])]
man = man.merge(cohort, on="key", how="inner").drop_duplicates("key")
man.to_csv(os.path.join(OUT, "manifesto_counts.csv"), index=False)
print("amostras do estudo com contagem:", len(man), "de", len(cohort))

def baixar(row):
    dest = os.path.join(OUT, row["key"] + ".tsv")
    if os.path.exists(dest) and os.path.getsize(dest) > 100000:
        return row["key"], "ok_existente", os.path.getsize(dest)
    u = "https://api.gdc.cancer.gov/data/" + row["file_id"]
    for tent in range(4):
        try:
            req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=300) as r, open(dest + ".part", "wb") as f:
                while True:
                    chunk = r.read(1024 * 256)
                    if not chunk:
                        break
                    f.write(chunk)
            os.replace(dest + ".part", dest)
            return row["key"], "ok", os.path.getsize(dest)
        except Exception as e:
            time.sleep(2 + tent * 3)
    return row["key"], "erro", 0

with ThreadPoolExecutor(max_workers=8) as ex:
    futs = [ex.submit(baixar, r) for _, r in man.iterrows()]
    for i, f in enumerate(as_completed(futs), 1):
        k, st, sz = f.result()
        if st != "ok_existente" or i % 25 == 0:
            print(i, k, st, sz, flush=True)
print("concluido")
