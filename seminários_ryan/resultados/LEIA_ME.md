# Resultados

## Figuras

| Grupo | Arquivos |
|---|---|
| Fenótipos | `figuras/01–06_*.png` |
| UMAP / heatmap / marcadores | `figuras/07–12_*.png` |
| ComBat | `figuras/13–14_*.png` |
| Fenótipos avançados | `figuras/15–25_*.png` |
| Grupos alternativos | `figuras/26–27_*.png` |
| Grupos separados | `figuras/28–35_*.png` |
| Volcano principais | `figuras/Volcano_01`, `Volcano_02`, `Volcano_03` (PNG) |
| Volcano complementares | `figuras/Volcano_04–12_*.png` |

> Todas as figuras são salvas apenas em **PNG**. Não há PDFs.

## Tabelas

| Grupo | Arquivos |
|---|---|
| Fenótipos e UMAP | `tabelas/01–12_*.csv` |
| Fenótipos avançados | `tabelas/13–20_*.csv` |
| Grupos alternativos | `tabelas/21–22_*.csv` |
| Grupos separados (recidivas, histologia, vírus, TSS) | `tabelas/23–35_*.csv` |
| Expressão diferencial completa | `tabelas/DE/DEG_*_completo.csv` |
| Apenas DEGs significativos | `tabelas/DE/DEG_*_significativos.csv` |
| Top 20 up/down | `tabelas/DE/TOP20_*.csv` |
| Top 5 genes por comparação | `tabelas/DE/RESUMO_TOP5_GENES.csv` |
| Resumo de DEGs | `tabelas/DE/RESUMO_DEGs_TODAS_COMPARACOES.csv` |
| ComBat antes/depois | `tabelas/DE/COMPARACAO_ComBat_antes_depois.csv` |
| Genes rotulados nos volcanos | `tabelas/ROTULOS_Volcano_*.csv` |

## Relatório consolidado

- `../relatorio_completo.md` — todos os docs + figuras em um único Markdown.

## Auditorias

- `../auditorias/AUDITORIA_DADOS.md`
- `../auditorias/AUDITORIA_VOLCANO_12.md`
- `../auditorias/BENCHMARKS.md`
- `../auditorias/COMPLIANCE.md`

## Controle de qualidade

| Arquivo | Conteúdo |
|---|---|
| `qc/log_00_preparar_dados.txt` | log da preparação |
| `qc/log_01_analise_exploratoria.txt` | log da exploratória |
| `qc/log_02_expressao_diferencial.txt` | log do limma |
| `qc/log_03_volcano.txt` | log dos volcanos |
| `qc/log_04_ComBat.txt` | log da demonstração do ComBat |
| `qc/log_05_fenotipos_avancados.txt` | log dos fenótipos avançados |
| `qc/log_06_grupos_alternativos.txt` | log dos grupos alternativos |
| `qc/resumo_grupos_separados.txt` | resumo das análises separadas |
| `qc/log_10_auditoria.txt` | log da auditoria |
| `qc/log_11_benchmarks.txt` | log dos benchmarks |
| `qc/log_12_compliance.txt` | log do compliance |
| `qc/resumo_*.txt` | resumos textuais |
| `qc/UMAP_ComBat_coordenadas.csv` | coordenadas UMAP antes/depois |

## Critério de significância

FDR (Benjamini-Hochberg) < 0,05 e |log2FC| > 1, salvo indicação contrária.
