# Índice dos resultados

Estrutura de saída do pipeline (`pipeline_hepato.R`) e dos scripts de análise.

```
results/
├── audit/         # Logs, benchmark, sessionInfo e relatórios de auditoria
├── deg/           # Genes diferencialmente expressos (limma)
├── enrichment/    # ORA (GO/KEGG/Reactome), GSEA e GSVA
├── figures/       # Heatmap dos top DEGs
├── ppi/           # Rede STRING (nós, arestas, topologia, hubs)
├── qc/            # Controle de qualidade (PCA, UMAP, sumário)
├── tables/        # Genes da via KEGG hsa05417
├── volcano/       # Volcano plot da via (PNG/PDF)
└── 3grupos/       # Análise em 3 grupos (Normal × Adjacente × LIHC)
```

## Descrição por pasta

| Pasta | Conteúdo principal |
|-------|--------------------|
| `audit/` | `pipeline_log.csv`, `benchmark_pipeline.csv`, `sessionInfo.txt`, relatórios `.md` |
| `deg/` | `DEG_LIHC_vs_Normal_full.csv`, `DEG_significant_LA.csv`, `top_up/down_genes.csv` |
| `enrichment/` | `GO_BP_*.csv`, `KEGG_*.csv`, `GSEA_*.csv`, `GSVA_*.csv`, `Reactome_*.csv` + figuras |
| `figures/` | `Heatmap_Top_DEGs.png` |
| `ppi/` | `PPI_STRING_nodes/edges_score700.csv`, `PPI_topology_metrics.csv`, `PPI_STRING_network.png` |
| `qc/` | `PCA_pre_analysis*.png`, `UMAP_pre_analysis.png`, `qc_summary.md` |
| `tables/` | `KEGG_hsa05417_genes.csv` (216 genes da via) |
| `volcano/` | `Volcano_LIHC_LA_pathway.png/pdf`, `Volcano_labeled_genes.csv` |
| `3grupos/` | `DEG_LIHC_vs_Adjacent.csv`, `DEG_Adjacent_vs_Normal.csv`, `KM_OS_*.png`, `PCA_3grupos*.png` |

## Principais figuras

- **Volcano plot:** `volcano/Volcano_LIHC_LA_pathway.png`
- **Rede PPI:** `ppi/PPI_STRING_network.png`
- **Heatmap:** `figures/Heatmap_Top_DEGs.png`
- **Enriquecimento:** `enrichment/enrichment_dotplot.png`, `enrichment/GSEA_KEGG_barplot.png`, `enrichment/Reactome_dotplot.png`
- **QC:** `qc/PCA_pre_analysis.png`, `qc/UMAP_pre_analysis.png`
- **3 grupos:** `3grupos/PCA_3grupos.png`
- **Kaplan-Meier:** `3grupos/KM_OS_*.png` (MMP1, CXCL2, MMP9, BAX, etc.)
