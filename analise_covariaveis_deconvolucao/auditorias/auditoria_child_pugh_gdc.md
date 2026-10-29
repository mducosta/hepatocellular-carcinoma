# Auditoria do Child-Pugh contra a fonte primária do GDC

Fonte: `dados\raw\gdc_clinico\cases_diagnoses.json`, consulta ao GDC em 28/10/2026.

## Resultado

- Tumores TCGA-LIHC no estudo: **369**.
- Código antigo (ausente tratado como B-C): **154** casos B-C (132 ausentes + 22 válidos).
- Código corrigido (ausente = ausência): **22** casos B-C.
- Xena sem valor utilizável: **132**; no GDC, **84** são `Unknown` e **48** estão vazios.
- Concordância entre Xena e GDC nos casos com valor conhecido: **100.0%** (n = 237).

## Crosstab Xena (linhas) × GDC (colunas)

| cp_xena   |   A |   B |   C |   Unknown |   ausente |
|:----------|----:|----:|----:|----------:|----------:|
| A         | 215 |   0 |   0 |         0 |         0 |
| B         |   0 |  21 |   0 |         0 |         0 |
| C         |   0 |   0 |   1 |         0 |         0 |
| ausente   |   0 |   0 |   0 |        84 |        48 |

## Interpretação

A codificação do GDC confirma a correção. O valor ausente na matriz Xena não é uma categoria clínica; parte corresponde a `Unknown` explícito (n = 84) e parte a campo vazio (n = 48). Tratar os dois como B-C produzia 154 casos B-C, contra 22 válidos. A concordância de 100% entre os valores conhecidos indica que a correção não alterou nenhum caso A, B ou C, apenas removeu a inflação por ausência.