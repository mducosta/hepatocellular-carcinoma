# Manifesto de renomeação

Em 26/10/2026, os arquivos antigos do seminário foram renomeados para nomes
descritivos. A pasta `dados/` não foi alterada. As análises avançadas e o artigo
final já tinham nomes descritivos e foram preservados. A pasta `artigos/`
(referências completas) é mantida apenas localmente e não é versionada.

## Regras aplicadas

| Antes | Depois | Exemplo |
|---|---|---|
| `README.md` | `LEIA_ME.md` | documentação principal |
| `relatorio_completo.md` | `relatorio_consolidado_seminario.md` | relatório consolidado |
| `requirements_R.txt` | `dependencias_R.txt` | dependências |
| `NN_descricao.md` (docs) | `descricao.md` | `02_ComBat.md` → `correcao_lote_ComBat.md` |
| `NN_descricao.R` (scripts) | `etapa_NN_descricao.R` | `02_expressao_diferencial.R` → `etapa_03_expressao_diferencial_limma.R` |
| `run_all.R` | `executar_todas_as_etapas.R` | executor do pipeline |
| `NN_descricao.png` | `figura_NN_descricao.png` | `01_amostras_por_grupo.png` → `figura_01_amostras_por_grupo.png` |
| `Volcano_NN...png` | `volcano_NN...png` | minúsculo e sem ambiguidade |
| `NN_descricao.csv` | `tabela_NN_descricao.csv` | `01_amostras_por_grupo.csv` → `tabela_01_amostras_por_grupo.csv` |
| `DEG_...` | `expressao_diferencial_...` | pasta `expressao_diferencial/` |
| `TOP20_...` | `top20_...` | top genes |
| `ROTULOS_Volcano_...` | `rotulos_volcano_...` | rótulos dos volcanos |
| `resultados/tabelas/DE/` | `resultados/tabelas/expressao_diferencial/` | pasta de DE |
| `resultados/qc/` | `resultados/controle_qualidade/` | logs e resumos |
| `AUDITORIA_DADOS.md` | `auditoria_dados.md` | auditoria |
| `AUDITORIA_VOLCANO_12.md` | `auditoria_volcano_ComBat.md` | auditoria do volcano 12 |
| `BENCHMARKS.md` | `benchmarks_pipeline.md` | benchmarks |
| `COMPLIANCE.md` | `compliance_pipeline.md` | compliance |

## Observações

- As referências internas em scripts, documentação, relatórios e tabelas foram
  atualizadas por um script local de renomeação.
- O script de renomeação não é versionado; este manifesto mantém o registro
  dos nomes antigos e novos.
- Arquivos que diferiam apenas por maiúsculas/minúsculas
  (`Volcano_` vs. `volcano_`, `ROTULOS_` vs. `rotulos_`,
  `AUDITORIA_DADOS.md` vs. `auditoria_dados.md`) exigiram renomeação em duas
  etapas no Windows, que é case-insensitive.
