# Auditoria de dados — TCGA-LIHC × GTEx

**Data da auditoria:** 2026-10-17

## 1. Resumo

|topico                                | resultado |detalhe                                                                         |
|:-------------------------------------|:---------:|:-------------------------------------------------------------------------------|
|Amostras na análise                   |    OK     |529 amostras; 369 tumor, 50 adjacente, 110 GTEx normal                          |
|Amostras totais (com recidivas)       |    OK     |531 amostras; recidivas excluídas: 2                                            |
|Duplicatas de amostra                 |    OK     |duplicatas: 0                                                                   |
|IDs TCGA/GTEx                         |    OK     |TCGA: 419; GTEx: 110                                                            |
|Amostras recidivadas (-02)            |    OK     |presentes no metadado completo: 2; excluídas da análise: 2                      |
|Matriz de expressão                   |    OK     |12401 genes x 529 amostras                                                      |
|Genes duplicados                      |    OK     |duplicatas: 0                                                                   |
|Valores ausentes na matriz            |    OK     |NAs: 0                                                                          |
|Valores infinitos na matriz           |    OK     |Inf: 0                                                                          |
|Amostras da matriz no metadado        |    OK     |interseção: 529                                                                 |
|Metadado completo na matriz           |    OK     |faltantes: 0                                                                    |
|Faixa de expressão                    |    OK     |min = -9.97; max = 17.3 (log2 TPM + 0,001)                                      |
|Sexo (TCGA tumor)                     |    OK     |faltantes: 0                                                                    |
|Idade (TCGA tumor)                    |    OK     |faltantes: 1; faixa: 16–90                                                      |
|Estágio (TCGA tumor)                  |    OK     |faltantes: 24; níveis: I, II, III, IV                                           |
|Grau (TCGA tumor)                     |    OK     |faltantes: 5; níveis: G1, G2, G3, G4                                            |
|Consistência estágio/estágio agrupado |    OK     |inconsistências: 0                                                              |
|Consistência grau/grau agrupado       |    OK     |inconsistências: 0                                                              |
|Sobrevida global (TCGA tumor)         |    OK     |OS/time disponíveis: 368; eventos: 128; tempos negativos: 0                     |
|OS = 1 com tempo > 0                  |    OK     |violações: 0                                                                    |
|Confundimento estudo x condição       |  ALERTA   |Cramér's V = 1; TCGA tumor e GTEx normal são perfeitamente separados por estudo |
|Arquivos essenciais                   |    OK     |ausentes: 0                                                                     |

## 2. Distribuição das amostras

|     | GTEx_Normal| TCGA_Adjacent| TCGA_Tumor|
|:----|-----------:|-------------:|----------:|
|GTEX |         110|             0|          0|
|TCGA |           0|            50|        369|

**Achado crítico:** a comparação TCGA × GTEx é estruturalmente confundida.
Todo tumor primário é TCGA e todo fígado normal é GTEx. Portanto, o efeito
de estudo/lote não é separável do efeito tumor vs. normal. Isso motivou:

- não aplicar ComBat na análise principal;
- usar a comparação interna TCGA tumor vs. adjacente como evidência primária;
- interpretar TCGA vs. GTEx com cautela.

## 3. Missingness clínico

|variavel          | n_na| n_total| percentual_na|
|:-----------------|----:|-------:|-------------:|
|vascular_invasion |  261|     369|          70.7|
|viral_status      |  206|     369|          55.8|
|fibrosis_ishak    |  159|     369|          43.1|
|fibrosis_advanced |  159|     369|          43.1|
|child_pugh        |  132|     369|          35.8|
|afp               |   93|     369|          25.2|
|afp_high          |   93|     369|          25.2|
|stage             |   24|     369|           6.5|
|ethnicity         |   19|     369|           5.1|
|race              |   10|     369|           2.7|
|grade             |    5|     369|           1.4|
|age               |    1|     369|           0.3|
|OS_time           |    1|     369|           0.3|
|sex               |    0|     369|           0.0|
|histology         |    0|     369|           0.0|
|OS                |    0|     369|           0.0|

## 4. Integridade da matriz

- Dimensões: 12401 genes × 529 amostras.
- NAs: 0; infinitos: 0.
- Duplicatas de gene: 0.
- Amostras da matriz ausentes no metadado: 0.

## 5. Conclusão

A base é consistente para as análises propostas, com três ressalvas:
confundimento estrutural TCGA × GTEx, missingness elevado em algumas
variáveis clínicas (invasão vascular, fibrose, AFP e sorologia viral) e
ausência de validação experimental. Nenhuma inconsistência crítica de
integridade foi detectada.
