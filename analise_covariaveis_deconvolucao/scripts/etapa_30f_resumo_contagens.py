#!/usr/bin/env python3
"""Consolida os resultados de DESeq2/edgeR (universo completo e comum) vs limma."""
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / "resultados" / "analises_avancadas" / "tabelas" / "sensibilidades"
DE = ROOT / "resultados" / "analises_avancadas" / "tabelas" / "expressao_diferencial_covariaveis"
FDR, LFC = 0.05, 1

lim_files = {
    "tumor_vs_adjacente": DE / "modelo_tumor_vs_adjacente_base_completo.csv",
    "estagio": DE / "modelo_estagio_base_completo.csv",
    "grau": DE / "modelo_grau_base_completo.csv",
}
lim = {k: pd.read_csv(v) for k, v in lim_files.items()}

rows = []
def add(univ, met, ct, f):
    d = pd.read_csv(f)
    d = d.dropna(subset=["gene"])
    ord_col = "baseMean" if "baseMean" in d.columns else "logCPM"
    d = d.sort_values(ord_col, ascending=False).drop_duplicates("gene")
    sig = d[d.regulation.isin(["Up", "Down"])]
    L = lim[ct]
    a = set(L.loc[L.regulation != "NS", "gene"]); b = set(sig.gene)
    com = set(L.gene) & set(d.gene)
    ov = a & b
    rows.append(dict(universo=univ, metodo=met, contraste=ct, n=int(len(d)),
                     genes=int(len(d)), degs=int(len(sig)), up=int((sig.regulation == "Up").sum()),
                     down=int((sig.regulation == "Down").sum()),
                     degs_limma=int(len(a)), degs_comuns=int(len(ov)),
                     jaccard=len(ov) / len(a | b) if a | b else np.nan,
                     spearman=pd.Series(L.set_index("gene").loc[list(com), "logFC"]).corr(
                         pd.Series(d.set_index("gene").loc[list(com), "log2FoldChange"]), method="spearman"),
                     direcao=float((np.sign(L.set_index("gene").loc[list(ov), "logFC"]).values ==
                                    np.sign(d.set_index("gene").loc[list(ov), "log2FoldChange"]).values).mean())
                     if ov else np.nan))

for met in ["DESeq2", "edgeR"]:
    for ct in lim_files:
        f = S / f"contagens_{met}_{ct}.csv"
        if f.exists():
            add("completo", met, ct, f)
for ct in lim_files:
    f = S / f"contagens_comum_DESeq2_{ct}.csv"
    if f.exists():
        add("comum_12401", "DESeq2", ct, f)

out = pd.DataFrame(rows)
out.to_csv(S / "contagens_resumo_consolidado.csv", index=False)
print(out.to_string(index=False))
