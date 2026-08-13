# Reprogramação Transcriptômica da Via de Lipídios e Aterosclerose no Hepatocarcinoma

Análise exploratória de **expressão diferencial** da via *Lipid and
Atherosclerosis* (KEGG `hsa05417`) entre **Hepatocarcinoma (LIHC, TCGA)** e
**fígado normal (GTEx)**, com integração de rede de interação proteína-proteína
(STRING), enriquecimento funcional (GO/KEGG) e controle de qualidade
pré-análise.

---

## Resumo

- **Desenho:** LIHC (TCGA) *vs.* Normal (GTEx Liver)
- **Via analisada:** KEGG `hsa05417` — Lipid and Atherosclerosis (216 genes; 212 na matriz)
- **Dados:** UCSC Xena (expressão gênica, escala log2)
- **Método de DE:** `limma` (FDR < 0,05 e |logFC| > 1)
- **Rede PPI:** STRING v12.0 (score ≥ 700)
- **Enriquecimento:** `clusterProfiler` (GO BP e KEGG)

### Principais resultados

| Classe | Contagem |
|--------|----------|
| Amostras de fígado | 481 (LIHC: 371 · Normal: 110) |
| Genes da via na matriz | 212 |
| Up-regulados em LIHC | 12 |
| Down-regulados em LIHC | 30 |
| Não significativos (NS) | 170 |

**Top DEGs (por significância):** CALML5, CXCL2, CALML6, IKBKB, TNFRSF10B,
**BAX**, HSP90AB1, STAT3, ABCA1, **CAMK2B**, **MMP1**, **MMP9**, CAMK2A
(mais detalhes em [`results/audit/RELATORIO_FINAL_LOCK_LIHC.md`](results/audit/RELATORIO_FINAL_LOCK_LIHC.md)).

---

## Estrutura do repositório

```
.
├── pipeline_hepato.R                 # Pipeline principal (reprodutível)
├── scripts/
│   └── análise_expressão_diferencial.R   # Versão original (legado)
├── data/
│   └── README.md                     # Instruções para obter os dados de entrada
├── results/
│   ├── audit/                        # Relatórios, logs, benchmark e sessionInfo
│   ├── deg/                          # Genes diferencialmente expressos (CSV)
│   ├── enrichment/                   # GO BP e KEGG (CSV)
│   ├── figures/                      # Heatmap dos top DEGs (PNG)
│   ├── ppi/                          # Rede STRING e métricas topológicas
│   ├── qc/                           # PCA, UMAP e sumário de QC pré-análise
│   ├── tables/                       # Genes da via hsa05417
│   └── volcano/                      # Volcano plot da via (PNG/PDF)
├── .gitignore
└── LICENSE
```

---

## Pré-requisitos

R ≥ 4.0 e os pacotes:

```
dplyr, tidyr, tibble, stringr, ggplot2, ggrepel, limma, edgeR, igraph,
ggraph, pheatmap, scales, clusterProfiler, org.Hs.eg.db, KEGGREST,
STRINGdb, readr, readxl, rio, uwot
```

O script instala automaticamente qualquer pacote ausente via `install.packages`
(CRAN) — os pacotes do Bioconductor (`limma`, `edgeR`, `clusterProfiler`,
`org.Hs.eg.db`, `STRINGdb`) são baixados pelo gerenciador do Bioconductor se
necessário.

---

## Como executar

1. **Baixe os dados** seguindo as instruções em [`data/README.md`](data/README.md)
   e salve como `data/liver_lip_aterosclerose.tsv`.

2. Execute o pipeline a partir da raiz do repositório:

   ```bash
   Rscript pipeline_hepato.R
   ```

   O script localiza a raiz do projeto automaticamente (inclusive a pasta `data/`),
   cria a estrutura `results/` e gera todos os outputs, o log
   (`results/audit/pipeline_log.csv`) e o benchmark
   (`results/audit/benchmark_pipeline.csv`).

---

## Reprodução dos resultados publicados

Os outputs já publicados em `results/` foram gerados com:

- **R** 4.6.0 (Windows 11 x64)
- Sessão completa registrada em [`results/audit/sessionInfo.txt`](results/audit/sessionInfo.txt)
- Tempo total de execução: ~462 s (ver [`results/audit/benchmark_summary.md`](results/audit/benchmark_summary.md))

---

## Documentação / auditoria

Relatórios de auditoria linha a linha, reanálise e verificação do resumo expandido
estão em [`results/audit/`](results/audit/):

- `RELATORIO_FINAL_LOCK_LIHC.md` — veredito final e resultados consolidados
- `AUDITORIA_CODIGO_LINHA_A_LINHA.md` — auditoria do código
- `AUDITORIA_REANALISE_LIHC.md` — auditoria da reanálise
- `AUDITORIA_RESUMO_EXPANDIDO.md` — cruzamento do resumo com código/dados

---

## Declaração de uso de inteligência artificial

> *Portaria CNPq nº 2.664/2026* — Este trabalho contou com o uso de ferramentas
> de inteligência artificial generativa (ChatGPT-5.5, OpenAI; DeepSeek-V4-Pro,
> DeepSeek) como suporte para processamento de dados transcriptômicos, revisão
> de código em linguagem R e editoração científica, sob supervisão e validação
> integral dos autores, que assumem total responsabilidade pelo conteúdo final.

---

## Limitações

- Estudo exploratório com dados secundários (TCGA/GTEx);
- Sem validação em coorte independente;
- Sem demonstração de causalidade;
- Sem validação de biomarcadores;
- Sem alegação de aplicabilidade clínica direta.

---

## Licença

[MIT](LICENSE) © 2026 Ryan de Paulo Santos.
