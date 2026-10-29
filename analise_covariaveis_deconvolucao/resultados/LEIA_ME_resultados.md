# Resultados

## Figuras

| Grupo | Arquivos |
|---|---|
| Fenótipos | `figuras/figura_01_*` a `figura_06_*` |
| UMAP, heatmap e marcadores | `figuras/figura_07_*` a `figura_12_*` |
| ComBat | `figuras/figura_13_*` e `figura_14_*` |
| Fenótipos avançados | `figuras/figura_15_*` a `figura_25_*` |
| Grupos alternativos | `figuras/figura_26_*` e `figura_27_*` |
| Grupos separados | `figuras/figura_28_*` a `figura_35_*` |
| Volcano principais | `figuras/volcano_01`, `volcano_02`, `volcano_03` |
| Volcano complementares | `figuras/volcano_04_*` a `volcano_12_*` |
| Análises avançadas | `analises_avancadas/figuras/` |
| Artigo final | `../artigo_final/figuras/` |

> As figuras de resultados são PNG em 300 dpi. O artigo final está em DOCX em
> `../artigo_final/`.

## Tabelas

| Grupo | Arquivos |
|---|---|
| Fenótipos e UMAP | `tabelas/tabela_01_*` a `tabela_12_*` |
| Fenótipos avançados | `tabelas/tabela_13_*` a `tabela_20_*` |
| Grupos alternativos | `tabelas/tabela_21_*` e `tabela_22_*` |
| Grupos separados | `tabelas/tabela_23_*` a `tabela_35_*` |
| Expressão diferencial completa | `tabelas/expressao_diferencial/expressao_diferencial_*_completo.csv` |
| Apenas DEGs significativos | `tabelas/expressao_diferencial/expressao_diferencial_*_significativos.csv` |
| Top 20 up/down | `tabelas/expressao_diferencial/top20_*.csv` |
| Top 5 genes por comparação | `tabelas/expressao_diferencial/resumo_top5_genes.csv` |
| Resumo de DEGs | `tabelas/expressao_diferencial/resumo_expressao_diferencial_todas_comparacoes.csv` |
| ComBat antes/depois | `tabelas/expressao_diferencial/comparacao_ComBat_antes_depois.csv` |
| Genes rotulados nos volcanos | `tabelas/rotulos_volcano_*.csv` |
| Análises avançadas | `analises_avancadas/tabelas/` |

## Relatório consolidado

- `../relatorio_consolidado_seminario.md` — documentos + figuras em Markdown.

## Auditorias

- `../auditorias/auditoria_analises_avancadas.md`
- `../auditorias/auditoria_dados.md`
- `../auditorias/auditoria_volcano_ComBat.md`
- `../auditorias/benchmarks_analises_avancadas.md`
- `../auditorias/compliance_analises_avancadas.md`
- `../auditorias/benchmarks_pipeline.md`
- `../auditorias/compliance_pipeline.md`

## Controle de qualidade

| Arquivo | Conteúdo |
|---|---|
| `controle_qualidade/log_etapa_00_download_dados.txt` | download |
| `controle_qualidade/log_etapa_01_preparar_matrizes_e_metadados.txt` | preparação |
| `controle_qualidade/log_etapa_02_analise_exploratoria.txt` | exploratória |
| `controle_qualidade/log_etapa_03_expressao_diferencial.txt` | limma |
| `controle_qualidade/log_etapa_04_volcano.txt` | volcanos |
| `controle_qualidade/log_etapa_05_ComBat.txt` | ComBat |
| `controle_qualidade/log_etapa_06_fenotipos_avancados.txt` | fenótipos avançados |
| `controle_qualidade/log_etapa_07_grupos_alternativos.txt` | grupos alternativos |
| `controle_qualidade/log_etapa_10_auditoria.txt` | auditoria |
| `controle_qualidade/log_etapa_11_benchmarks.txt` | benchmarks |
| `controle_qualidade/log_etapa_12_compliance.txt` | compliance |
| `controle_qualidade/resumo_*.txt` | resumos textuais |
| `controle_qualidade/coordenadas_UMAP_ComBat.csv` | coordenadas UMAP |

## Critério de significância

FDR (Benjamini-Hochberg) < 0,05 e |log2FC| > 1, salvo indicação contrária.
