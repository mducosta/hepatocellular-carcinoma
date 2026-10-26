# Auditoria das análises avançadas — TCGA-LIHC × GTEx

Gerado em 2026-10-26 14:26.

## 1. Inventário e integridade

| item                                          | status   | evidencia                                                                               |
|:----------------------------------------------|:---------|:----------------------------------------------------------------------------------------|
| Total de amostras no metadado                 | OK       | 529 amostras                                                                            |
| Tumor / adjacente / GTEx                      | OK       | 369 / 50 / 110                                                                          |
| Tumores recorrentes                           | INFO     | 2 recidivas foram excluídas na preparação; 0 no metadado analítico                      |
| GTEx com mesmo doador em >1 amostra           | ATENÇÃO  | 1 doador(es); 2 amostras; análise principal usa 1 por doador                            |
| Pacientes com tumor e adjacente               | OK       | 50 pares                                                                                |
| Child-Pugh ausente                            | ATENÇÃO  | 132 tumores sem Child-Pugh; a versão anterior recodificava ausentes como B-C            |
| Sorologia desconhecida                        | ATENÇÃO  | 206 tumores; não é 'sem vírus'                                                          |
| Histologia 'Other'                            | ATENÇÃO  | 10 tumores (7 mistos + 3 fibrolamelares)                                                |
| Estágio IV                                    | ATENÇÃO  | 5 tumores                                                                               |
| Grau G4                                       | ATENÇÃO  | 11 tumores                                                                              |
| Idade GTEx                                    | ATENÇÃO  | faixa etária convertida ao ponto médio                                                  |
| Dados de expressão                            | OK       | log2(TPM + 0,001); 12.401 genes; 528 amostras após deduplicação                         |
| Matriz de expressão                           | OK       | 12401 genes × 529 amostras; análise principal com 528 após deduplicação                 |
| Valores ausentes na matriz                    | OK       | 0 valores ausentes                                                                      |
| Genes com variância zero                      | OK       | 0 genes                                                                                 |
| Escala log2                                   | OK       | min=-9.966; max=17.304; mediana=3.072                                                   |
| EPIC: convergência                            | ATENÇÃO  | 84 de 528 amostras com código != 0                                                      |
| EPIC otherCells vs marcadores de hepatócito   | ATENÇÃO  | Spearman negativo com ALB/TTR; otherCells NÃO deve ser chamado de fração de hepatócitos |
| xCell Hepatocytes vs marcadores de hepatócito | OK       | Spearman positivo com ALB, HP, TTR e CYP3A4                                             |
| Concordância entre métodos                    | ATENÇÃO  | Spearman mediano entre métodos = 0.01                                                   |
| Colinearidade dos escores de composição       | INFO     | xCell_Hepatocytes: VIF=1.36; xCell_ImmuneScore: VIF=1.05; xCell_StromaScore: VIF=1.38   |

## 2. Associações e confundimento

| var_a                      | var_b                      |   n |   cramers_v |     p |
|:---------------------------|:---------------------------|----:|------------:|------:|
| group                      | study                      | 529 |       1.000 | 0.000 |
| stage_group2               | tss                        | 345 |       0.520 | 0.000 |
| grade_group2               | tss                        | 364 |       0.430 | 0.002 |
| histology                  | tss                        | 369 |       0.414 | 0.000 |
| child_pugh_group_corrigido | tss                        | 237 |       0.324 | 0.408 |
| sex                        | tss                        | 369 |       0.321 | 0.424 |
| grade_group2               | afp_group                  | 364 |       0.320 | 0.000 |
| stage_group2               | afp_group                  | 345 |       0.284 | 0.000 |
| stage_group2               | viral_status               | 345 |       0.276 | 0.000 |
| grade_group2               | viral_status               | 364 |       0.198 | 0.003 |
| stage_group2               | child_pugh_group_corrigido | 227 |       0.101 | 0.127 |
| sex                        | group                      | 529 |       0.078 | 0.200 |
| sex                        | stage_group2               | 345 |       0.075 | 0.163 |
| sex                        | grade_group2               | 364 |       0.044 | 0.404 |
| grade_group2               | child_pugh_group_corrigido | 236 |       0.016 | 0.801 |

### TSS com maior proporção de estágio III–IV

| tss   |   I-II |   III-IV |   n |   pct_III_IV |   child_pugh_ausente |
|:------|-------:|---------:|----:|-------------:|---------------------:|
| DD    |    119 |       25 | 144 |         17.4 |                    8 |
| CC    |      6 |       26 |  32 |         81.2 |                   32 |
| G3    |     26 |        6 |  32 |         18.8 |                    1 |
| 2Y    |     17 |        1 |  18 |          5.6 |                   14 |
| ED    |      8 |        7 |  15 |         46.7 |                   15 |
| BC    |      6 |        6 |  12 |         50.0 |                   19 |
| FV    |     10 |        1 |  11 |          9.1 |                   11 |
| RC    |      7 |        2 |   9 |         22.2 |                    1 |
| UB    |      6 |        2 |   8 |         25.0 |                    0 |
| EP    |      5 |        2 |   7 |         28.6 |                    2 |
| XR    |      3 |        2 |   5 |         40.0 |                    2 |
| K7    |      4 |        0 |   4 |          0.0 |                    3 |

## 3. Missingness corrigida (tumor)

| variavel          |   ausente |   total |   percentual |
|:------------------|----------:|--------:|-------------:|
| vascular_invasion |       261 |     369 |       70.700 |
| viral_status      |       206 |     369 |       55.800 |
| fibrosis_ishak    |       159 |     369 |       43.100 |
| child_pugh        |       132 |     369 |       35.800 |
| afp               |        93 |     369 |       25.200 |
| stage             |        24 |     369 |        6.500 |
| grade             |         5 |     369 |        1.400 |
| OS_time           |         1 |     369 |        0.300 |
| age               |         1 |     369 |        0.300 |
| sex               |         0 |     369 |        0.000 |
| histology         |         0 |     369 |        0.000 |
| tss               |         0 |     369 |        0.000 |
| OS                |         0 |     369 |        0.000 |

## 4. Fluxo de casos

| modelo                                                  |   n_total |   n_tumor |   n_adjacente |
|:--------------------------------------------------------|----------:|----------:|--------------:|
| Tumor vs adjacente: base (sexo, idade)                  |       418 |       368 |        50.000 |
| Tumor vs adjacente: + estágio e grau                    |       382 |       342 |        40.000 |
| Tumor vs adjacente: + sorologia, Child-Pugh, histologia |       257 |       226 |        31.000 |
| Tumor vs adjacente: + TSS                               |       257 |       226 |        31.000 |
| Estágio: base                                           |       344 |       344 |       nan     |
| Estágio: + grau e clínica                               |       226 |       226 |       nan     |
| Estágio: + TSS                                          |       226 |       226 |       nan     |
| Grau: base                                              |       363 |       363 |       nan     |
| Grau: + estágio e clínica                               |       226 |       226 |       nan     |
| Grau: + TSS                                             |       226 |       226 |       nan     |

## 5. Atenuação de sinais

| comparacao                                         |   degs_base |   degs_alt |   jaccard |   atenuacao_mediana |   n_base |   n_alt |
|:---------------------------------------------------|------------:|-----------:|----------:|--------------------:|---------:|--------:|
| Estágio: coorte completa vs mesmos casos completos |          46 |          0 |     0.000 |               0.650 |    12401 |   12401 |
| Estágio: sem vs com TSS                            |          46 |          0 |     0.000 |               0.467 |    12401 |   12401 |
| Estágio: excluindo CC/ED                           |          46 |          0 |     0.000 |               0.408 |    12401 |   12401 |
| Grau: base vs clínico + TSS                        |         172 |        247 |     0.387 |              -0.028 |    12401 |   12401 |
| Tumor vs adjacente: não pareado vs pareado         |        4130 |       3710 |     0.764 |               0.052 |    12401 |   12401 |
| Tumor vs adjacente: sem vs com composição          |        4130 |       1897 |     0.317 |               0.368 |    12401 |   12401 |
| Tumor vs GTEx: sem vs com composição               |        3391 |       2732 |     0.554 |               0.187 |    12401 |   12401 |
| Grau: sem vs com composição                        |         172 |          0 |     0.000 |               0.692 |    12401 |   12401 |
| Estágio: sem vs com composição                     |          46 |          0 |     0.000 |               0.597 |    12401 |   12401 |

## 6. QC da deconvolução

| celula                  | marcador   |   spearman | familia        |
|:------------------------|:-----------|-----------:|:---------------|
| EPIC_otherCells         | ALB        |     -0.327 | Hepatocytes    |
| EPIC_otherCells         | APOA1      |     -0.379 | Hepatocytes    |
| EPIC_otherCells         | HP         |     -0.153 | Hepatocytes    |
| EPIC_otherCells         | TTR        |     -0.509 | Hepatocytes    |
| EPIC_otherCells         | CYP3A4     |     -0.306 | Hepatocytes    |
| EPIC_otherCells         | CYP2E1     |     -0.084 | Hepatocytes    |
| EPIC_otherCells         | KRT19      |      0.290 | Cholangiocytes |
| EPIC_otherCells         | EPCAM      |      0.262 | Cholangiocytes |
| EPIC_otherCells         | KRT7       |      0.042 | Cholangiocytes |
| EPIC_otherCells         | SOX9       |      0.139 | Cholangiocytes |
| EPIC_otherCells         | PECAM1     |     -0.012 | Endothelial    |
| EPIC_otherCells         | CDH5       |     -0.006 | Endothelial    |
| EPIC_otherCells         | VWF        |     -0.180 | Endothelial    |
| EPIC_otherCells         | KDR        |     -0.014 | Endothelial    |
| EPIC_otherCells         | COL1A1     |      0.173 | Fibroblasts    |
| EPIC_otherCells         | COL1A2     |      0.145 | Fibroblasts    |
| EPIC_otherCells         | ACTA2      |     -0.019 | Fibroblasts    |
| EPIC_otherCells         | PDGFRB     |      0.133 | Fibroblasts    |
| EPIC_otherCells         | CD68       |      0.192 | Macrophages    |
| EPIC_otherCells         | CD14       |     -0.372 | Macrophages    |
| EPIC_otherCells         | CSF1R      |      0.178 | Macrophages    |
| EPIC_otherCells         | AIF1       |      0.016 | Macrophages    |
| EPIC_otherCells         | CD3D       |     -0.040 | T cells        |
| EPIC_otherCells         | CD3E       |     -0.019 | T cells        |
| EPIC_otherCells         | CD2        |     -0.035 | T cells        |
| EPIC_otherCells         | TRAC       |     -0.025 | T cells        |
| EPIC_otherCells         | CD79A      |      0.078 | B cells        |
| EPIC_otherCells         | NKG7       |     -0.113 | NK             |
| EPIC_otherCells         | GNLY       |      0.057 | NK             |
| EPIC_otherCells         | S100A8     |      0.293 | Neutrophils    |
| EPIC_otherCells         | S100A9     |      0.397 | Neutrophils    |
| xCell_Hepatocytes       | ALB        |      0.474 | Hepatocytes    |
| xCell_Hepatocytes       | APOA1      |      0.448 | Hepatocytes    |
| xCell_Hepatocytes       | HP         |      0.623 | Hepatocytes    |
| xCell_Hepatocytes       | TTR        |      0.616 | Hepatocytes    |
| xCell_Hepatocytes       | CYP3A4     |      0.643 | Hepatocytes    |
| xCell_Hepatocytes       | CYP2E1     |      0.415 | Hepatocytes    |
| xCell_Hepatocytes       | KRT19      |     -0.085 | Cholangiocytes |
| xCell_Hepatocytes       | EPCAM      |     -0.206 | Cholangiocytes |
| xCell_Hepatocytes       | KRT7       |      0.008 | Cholangiocytes |
| xCell_Hepatocytes       | SOX9       |     -0.483 | Cholangiocytes |
| xCell_Hepatocytes       | PECAM1     |     -0.414 | Endothelial    |
| xCell_Hepatocytes       | CDH5       |     -0.146 | Endothelial    |
| xCell_Hepatocytes       | VWF        |     -0.342 | Endothelial    |
| xCell_Hepatocytes       | KDR        |      0.067 | Endothelial    |
| xCell_Hepatocytes       | COL1A1     |     -0.380 | Fibroblasts    |
| xCell_Hepatocytes       | COL1A2     |     -0.374 | Fibroblasts    |
| xCell_Hepatocytes       | ACTA2      |     -0.273 | Fibroblasts    |
| xCell_Hepatocytes       | PDGFRB     |     -0.459 | Fibroblasts    |
| xCell_Hepatocytes       | CD68       |     -0.329 | Macrophages    |
| xCell_Hepatocytes       | CD14       |      0.460 | Macrophages    |
| xCell_Hepatocytes       | CSF1R      |     -0.187 | Macrophages    |
| xCell_Hepatocytes       | AIF1       |     -0.327 | Macrophages    |
| xCell_Hepatocytes       | CD3D       |     -0.448 | T cells        |
| xCell_Hepatocytes       | CD3E       |     -0.224 | T cells        |
| xCell_Hepatocytes       | CD2        |     -0.314 | T cells        |
| xCell_Hepatocytes       | TRAC       |     -0.334 | T cells        |
| xCell_Hepatocytes       | CD79A      |     -0.050 | B cells        |
| xCell_Hepatocytes       | NKG7       |      0.028 | NK             |
| xCell_Hepatocytes       | GNLY       |      0.036 | NK             |
| xCell_Hepatocytes       | S100A8     |      0.241 | Neutrophils    |
| xCell_Hepatocytes       | S100A9     |      0.040 | Neutrophils    |
| xCell_Endothelial cells | ALB        |      0.067 | Hepatocytes    |
| xCell_Endothelial cells | APOA1      |      0.011 | Hepatocytes    |
| xCell_Endothelial cells | HP         |     -0.021 | Hepatocytes    |
| xCell_Endothelial cells | TTR        |      0.039 | Hepatocytes    |
| xCell_Endothelial cells | CYP3A4     |     -0.059 | Hepatocytes    |
| xCell_Endothelial cells | CYP2E1     |     -0.097 | Hepatocytes    |
| xCell_Endothelial cells | KRT19      |     -0.105 | Cholangiocytes |
| xCell_Endothelial cells | EPCAM      |     -0.031 | Cholangiocytes |
| xCell_Endothelial cells | KRT7       |      0.054 | Cholangiocytes |
| xCell_Endothelial cells | SOX9       |      0.163 | Cholangiocytes |
| xCell_Endothelial cells | PECAM1     |      0.754 | Endothelial    |
| xCell_Endothelial cells | CDH5       |      0.705 | Endothelial    |
| xCell_Endothelial cells | VWF        |      0.766 | Endothelial    |
| xCell_Endothelial cells | KDR        |      0.586 | Endothelial    |
| xCell_Endothelial cells | COL1A1     |      0.402 | Fibroblasts    |
| xCell_Endothelial cells | COL1A2     |      0.445 | Fibroblasts    |
| xCell_Endothelial cells | ACTA2      |      0.508 | Fibroblasts    |
| xCell_Endothelial cells | PDGFRB     |      0.439 | Fibroblasts    |
| xCell_Endothelial cells | CD68       |      0.331 | Macrophages    |
| xCell_Endothelial cells | CD14       |      0.270 | Macrophages    |
| xCell_Endothelial cells | CSF1R      |      0.324 | Macrophages    |
| xCell_Endothelial cells | AIF1       |      0.418 | Macrophages    |
| xCell_Endothelial cells | CD3D       |      0.310 | T cells        |
| xCell_Endothelial cells | CD3E       |      0.305 | T cells        |
| xCell_Endothelial cells | CD2        |      0.342 | T cells        |
| xCell_Endothelial cells | TRAC       |      0.359 | T cells        |
| xCell_Endothelial cells | CD79A      |      0.131 | B cells        |
| xCell_Endothelial cells | NKG7       |      0.083 | NK             |
| xCell_Endothelial cells | GNLY       |     -0.068 | NK             |
| xCell_Endothelial cells | S100A8     |     -0.341 | Neutrophils    |
| xCell_Endothelial cells | S100A9     |     -0.296 | Neutrophils    |
| EPIC_CAFs               | ALB        |     -0.041 | Hepatocytes    |
| EPIC_CAFs               | APOA1      |     -0.156 | Hepatocytes    |
| EPIC_CAFs               | HP         |     -0.161 | Hepatocytes    |
| EPIC_CAFs               | TTR        |     -0.235 | Hepatocytes    |
| EPIC_CAFs               | CYP3A4     |     -0.171 | Hepatocytes    |
| EPIC_CAFs               | CYP2E1     |     -0.226 | Hepatocytes    |
| EPIC_CAFs               | KRT19      |      0.458 | Cholangiocytes |
| EPIC_CAFs               | EPCAM      |      0.410 | Cholangiocytes |
| EPIC_CAFs               | KRT7       |      0.362 | Cholangiocytes |
| EPIC_CAFs               | SOX9       |      0.278 | Cholangiocytes |
| EPIC_CAFs               | PECAM1     |      0.431 | Endothelial    |
| EPIC_CAFs               | CDH5       |      0.437 | Endothelial    |
| EPIC_CAFs               | VWF        |      0.112 | Endothelial    |
| EPIC_CAFs               | KDR        |      0.305 | Endothelial    |
| EPIC_CAFs               | COL1A1     |      0.728 | Fibroblasts    |
| EPIC_CAFs               | COL1A2     |      0.721 | Fibroblasts    |
| EPIC_CAFs               | ACTA2      |      0.509 | Fibroblasts    |
| EPIC_CAFs               | PDGFRB     |      0.533 | Fibroblasts    |
| EPIC_CAFs               | CD68       |      0.340 | Macrophages    |
| EPIC_CAFs               | CD14       |     -0.084 | Macrophages    |
| EPIC_CAFs               | CSF1R      |      0.423 | Macrophages    |
| EPIC_CAFs               | AIF1       |      0.320 | Macrophages    |
| EPIC_CAFs               | CD3D       |      0.266 | T cells        |
| EPIC_CAFs               | CD3E       |      0.355 | T cells        |
| EPIC_CAFs               | CD2        |      0.341 | T cells        |
| EPIC_CAFs               | TRAC       |      0.325 | T cells        |
| EPIC_CAFs               | CD79A      |      0.319 | B cells        |
| EPIC_CAFs               | NKG7       |      0.099 | NK             |
| EPIC_CAFs               | GNLY       |      0.133 | NK             |
| EPIC_CAFs               | S100A8     |      0.179 | Neutrophils    |
| EPIC_CAFs               | S100A9     |      0.256 | Neutrophils    |
| EPIC_Macrophages        | ALB        |      0.372 | Hepatocytes    |
| EPIC_Macrophages        | APOA1      |      0.424 | Hepatocytes    |
| EPIC_Macrophages        | HP         |      0.204 | Hepatocytes    |
| EPIC_Macrophages        | TTR        |      0.558 | Hepatocytes    |
| EPIC_Macrophages        | CYP3A4     |      0.340 | Hepatocytes    |
| EPIC_Macrophages        | CYP2E1     |      0.126 | Hepatocytes    |
| EPIC_Macrophages        | KRT19      |     -0.323 | Cholangiocytes |
| EPIC_Macrophages        | EPCAM      |     -0.275 | Cholangiocytes |
| EPIC_Macrophages        | KRT7       |     -0.093 | Cholangiocytes |
| EPIC_Macrophages        | SOX9       |     -0.161 | Cholangiocytes |
| EPIC_Macrophages        | PECAM1     |     -0.072 | Endothelial    |
| EPIC_Macrophages        | CDH5       |     -0.079 | Endothelial    |
| EPIC_Macrophages        | VWF        |      0.133 | Endothelial    |
| EPIC_Macrophages        | KDR        |     -0.054 | Endothelial    |
| EPIC_Macrophages        | COL1A1     |     -0.248 | Fibroblasts    |
| EPIC_Macrophages        | COL1A2     |     -0.223 | Fibroblasts    |
| EPIC_Macrophages        | ACTA2      |     -0.070 | Fibroblasts    |
| EPIC_Macrophages        | PDGFRB     |     -0.198 | Fibroblasts    |
| EPIC_Macrophages        | CD68       |     -0.214 | Macrophages    |
| EPIC_Macrophages        | CD14       |      0.380 | Macrophages    |
| EPIC_Macrophages        | CSF1R      |     -0.219 | Macrophages    |
| EPIC_Macrophages        | AIF1       |     -0.057 | Macrophages    |
| EPIC_Macrophages        | CD3D       |      0.003 | T cells        |
| EPIC_Macrophages        | CD3E       |     -0.039 | T cells        |
| EPIC_Macrophages        | CD2        |     -0.024 | T cells        |
| EPIC_Macrophages        | TRAC       |     -0.026 | T cells        |
| EPIC_Macrophages        | CD79A      |     -0.116 | B cells        |
| EPIC_Macrophages        | NKG7       |      0.077 | NK             |
| EPIC_Macrophages        | GNLY       |     -0.079 | NK             |
| EPIC_Macrophages        | S100A8     |     -0.293 | Neutrophils    |
| EPIC_Macrophages        | S100A9     |     -0.394 | Neutrophils    |
| xCell_Macrophages       | ALB        |      0.195 | Hepatocytes    |
| xCell_Macrophages       | APOA1      |      0.143 | Hepatocytes    |
| xCell_Macrophages       | HP         |      0.538 | Hepatocytes    |
| xCell_Macrophages       | TTR        |      0.173 | Hepatocytes    |
| xCell_Macrophages       | CYP3A4     |      0.152 | Hepatocytes    |
| xCell_Macrophages       | CYP2E1     |      0.109 | Hepatocytes    |
| xCell_Macrophages       | KRT19      |      0.142 | Cholangiocytes |
| xCell_Macrophages       | EPCAM      |      0.049 | Cholangiocytes |
| xCell_Macrophages       | KRT7       |      0.152 | Cholangiocytes |
| xCell_Macrophages       | SOX9       |     -0.153 | Cholangiocytes |
| xCell_Macrophages       | PECAM1     |      0.081 | Endothelial    |
| xCell_Macrophages       | CDH5       |      0.168 | Endothelial    |
| xCell_Macrophages       | VWF        |      0.016 | Endothelial    |
| xCell_Macrophages       | KDR        |      0.221 | Endothelial    |
| xCell_Macrophages       | COL1A1     |      0.137 | Fibroblasts    |
| xCell_Macrophages       | COL1A2     |      0.158 | Fibroblasts    |
| xCell_Macrophages       | ACTA2      |      0.062 | Fibroblasts    |
| xCell_Macrophages       | PDGFRB     |     -0.120 | Fibroblasts    |
| xCell_Macrophages       | CD68       |      0.370 | Macrophages    |
| xCell_Macrophages       | CD14       |      0.430 | Macrophages    |
| xCell_Macrophages       | CSF1R      |      0.476 | Macrophages    |
| xCell_Macrophages       | AIF1       |      0.384 | Macrophages    |
| xCell_Macrophages       | CD3D       |      0.090 | T cells        |
| xCell_Macrophages       | CD3E       |      0.259 | T cells        |
| xCell_Macrophages       | CD2        |      0.200 | T cells        |
| xCell_Macrophages       | TRAC       |      0.178 | T cells        |
| xCell_Macrophages       | CD79A      |      0.294 | B cells        |
| xCell_Macrophages       | NKG7       |      0.210 | NK             |
| xCell_Macrophages       | GNLY       |      0.051 | NK             |
| xCell_Macrophages       | S100A8     |      0.283 | Neutrophils    |
| xCell_Macrophages       | S100A9     |      0.289 | Neutrophils    |
| EPIC_CD8_Tcells         | ALB        |      0.089 | Hepatocytes    |
| EPIC_CD8_Tcells         | APOA1      |      0.114 | Hepatocytes    |
| EPIC_CD8_Tcells         | HP         |     -0.030 | Hepatocytes    |
| EPIC_CD8_Tcells         | TTR        |      0.148 | Hepatocytes    |
| EPIC_CD8_Tcells         | CYP3A4     |      0.107 | Hepatocytes    |
| EPIC_CD8_Tcells         | CYP2E1     |     -0.038 | Hepatocytes    |
| EPIC_CD8_Tcells         | KRT19      |     -0.078 | Cholangiocytes |
| EPIC_CD8_Tcells         | EPCAM      |     -0.079 | Cholangiocytes |
| EPIC_CD8_Tcells         | KRT7       |     -0.053 | Cholangiocytes |
| EPIC_CD8_Tcells         | SOX9       |      0.037 | Cholangiocytes |
| EPIC_CD8_Tcells         | PECAM1     |      0.016 | Endothelial    |
| EPIC_CD8_Tcells         | CDH5       |     -0.009 | Endothelial    |
| EPIC_CD8_Tcells         | VWF        |      0.046 | Endothelial    |
| EPIC_CD8_Tcells         | KDR        |      0.009 | Endothelial    |
| EPIC_CD8_Tcells         | COL1A1     |     -0.091 | Fibroblasts    |
| EPIC_CD8_Tcells         | COL1A2     |     -0.091 | Fibroblasts    |
| EPIC_CD8_Tcells         | ACTA2      |     -0.001 | Fibroblasts    |
| EPIC_CD8_Tcells         | PDGFRB     |     -0.085 | Fibroblasts    |
| EPIC_CD8_Tcells         | CD68       |     -0.032 | Macrophages    |
| EPIC_CD8_Tcells         | CD14       |      0.053 | Macrophages    |
| EPIC_CD8_Tcells         | CSF1R      |     -0.059 | Macrophages    |
| EPIC_CD8_Tcells         | AIF1       |      0.008 | Macrophages    |
| EPIC_CD8_Tcells         | CD3D       |     -0.002 | T cells        |
| EPIC_CD8_Tcells         | CD3E       |     -0.036 | T cells        |
| EPIC_CD8_Tcells         | CD2        |     -0.020 | T cells        |
| EPIC_CD8_Tcells         | TRAC       |     -0.021 | T cells        |
| EPIC_CD8_Tcells         | CD79A      |     -0.132 | B cells        |
| EPIC_CD8_Tcells         | NKG7       |      0.018 | NK             |
| EPIC_CD8_Tcells         | GNLY       |     -0.060 | NK             |
| EPIC_CD8_Tcells         | S100A8     |     -0.178 | Neutrophils    |
| EPIC_CD8_Tcells         | S100A9     |     -0.226 | Neutrophils    |
| xCell_CD8+ T-cells      | ALB        |      0.014 | Hepatocytes    |
| xCell_CD8+ T-cells      | APOA1      |     -0.067 | Hepatocytes    |
| xCell_CD8+ T-cells      | HP         |     -0.081 | Hepatocytes    |
| xCell_CD8+ T-cells      | TTR        |     -0.074 | Hepatocytes    |
| xCell_CD8+ T-cells      | CYP3A4     |     -0.063 | Hepatocytes    |
| xCell_CD8+ T-cells      | CYP2E1     |     -0.058 | Hepatocytes    |
| xCell_CD8+ T-cells      | KRT19      |      0.078 | Cholangiocytes |
| xCell_CD8+ T-cells      | EPCAM      |      0.004 | Cholangiocytes |
| xCell_CD8+ T-cells      | KRT7       |      0.058 | Cholangiocytes |
| xCell_CD8+ T-cells      | SOX9       |      0.098 | Cholangiocytes |
| xCell_CD8+ T-cells      | PECAM1     |      0.098 | Endothelial    |
| xCell_CD8+ T-cells      | CDH5       |      0.058 | Endothelial    |
| xCell_CD8+ T-cells      | VWF        |      0.065 | Endothelial    |
| xCell_CD8+ T-cells      | KDR        |      0.049 | Endothelial    |
| xCell_CD8+ T-cells      | COL1A1     |      0.098 | Fibroblasts    |
| xCell_CD8+ T-cells      | COL1A2     |      0.119 | Fibroblasts    |
| xCell_CD8+ T-cells      | ACTA2      |      0.066 | Fibroblasts    |
| xCell_CD8+ T-cells      | PDGFRB     |      0.050 | Fibroblasts    |
| xCell_CD8+ T-cells      | CD68       |      0.095 | Macrophages    |
| xCell_CD8+ T-cells      | CD14       |      0.035 | Macrophages    |
| xCell_CD8+ T-cells      | CSF1R      |      0.167 | Macrophages    |
| xCell_CD8+ T-cells      | AIF1       |      0.172 | Macrophages    |
| xCell_CD8+ T-cells      | CD3D       |      0.142 | T cells        |
| xCell_CD8+ T-cells      | CD3E       |      0.160 | T cells        |
| xCell_CD8+ T-cells      | CD2        |      0.168 | T cells        |
| xCell_CD8+ T-cells      | TRAC       |      0.150 | T cells        |
| xCell_CD8+ T-cells      | CD79A      |      0.177 | B cells        |
| xCell_CD8+ T-cells      | NKG7       |      0.195 | NK             |
| xCell_CD8+ T-cells      | GNLY       |      0.080 | NK             |
| xCell_CD8+ T-cells      | S100A8     |      0.039 | Neutrophils    |
| xCell_CD8+ T-cells      | S100A9     |      0.067 | Neutrophils    |

## 7. Concordância entre métodos

| conceito          | metodo_a         | metodo_b                | metodo_c                 |   spearman_ab |   spearman_ac |   spearman_bc |
|:------------------|:-----------------|:------------------------|:-------------------------|--------------:|--------------:|--------------:|
| Hepatócitos/other | EPIC_otherCells  | xCell_Hepatocytes       | nan                      |        -0.306 |       nan     |       nan     |
| Endotélio         | EPIC_Endothelial | xCell_Endothelial cells | nan                      |         0.227 |       nan     |       nan     |
| Fibroblastos/CAFs | EPIC_CAFs        | xCell_Fibroblasts       | nan                      |         0.156 |       nan     |       nan     |
| Macrófagos        | EPIC_Macrophages | xCell_Macrophages       | quanTIseq_Macrophages.M1 |         0.105 |        -0.397 |         0.254 |
| CD8 T             | EPIC_CD8_Tcells  | xCell_CD8+ T-cells      | quanTIseq_T.cells.CD8    |        -0.018 |        -0.064 |         0.216 |
| CD4 T             | EPIC_CD4_Tcells  | xCell_CD4+ T-cells      | quanTIseq_T.cells.CD4    |        -0.026 |        -0.102 |        -0.023 |
| B                 | EPIC_Bcells      | xCell_B-cells           | quanTIseq_B.cells        |         0.110 |         0.081 |         0.207 |
| NK                | EPIC_NKcells     | xCell_NK cells          | quanTIseq_NK.cells       |         0.012 |        -0.055 |         0.017 |

## 8. Achados críticos

1. **Confundimento grupo × estudo.** Cramér's V = 1,0; logFC tumor vs GTEx não separa biologia de lote.
2. **Child-Pugh mal codificado.** 132 tumores sem Child-Pugh foram tratados como B-C na versão anterior; aqui são NA. O n de B-C cai de 154 para 22.
3. **Estágio confundido com TSS.** CC concentra 26/90 tumores III-IV e 100% dos CC sem Child-Pugh; o ajuste por TSS zera os DEGs de estágio.
4. **Estágio frágil à restrição amostral.** 46 DEGs na coorte completa, 0 nos mesmos casos completos; logFC mediano cai de 1,31 para 0,52.
5. **Grau robusto, mas composicional.** 172–337 DEGs em modelos diferentes; ajuste por composição zera os DEGs.
6. **Tumor vs adjacente robusto.** 4.130 DEGs no basal; 3.710 no pareado e 3.740 no clínico; Jaccard 0,76–0,83.
7. **Composição explica parte do sinal.** Tumor vs adjacente cai de 4.130 para 1.897 DEGs após ajuste por xCell.
8. **Sexo não é só composição.** 13 DEGs no GTEx e 175 no tumor; ajuste por composição reduz tumor para 90, mas não zera.
9. **Histologia não interpretável.** Apenas 10 casos 'Other'; 740–1.020 DEGs são instáveis e não devem ser interpretados como biologia.
10. **Sorologia/Child-Pugh/fibrose/invasão.** Zero DEGs com FDR < 0,05 e |logFC| > 1; poder baixo e missingness alta.
11. **Deconvolução não é fração verdadeira.** EPIC otherCells não é hepatócito; xCell é escore de enriquecimento; concordância entre métodos é baixa.

## 9. Conclusão da auditoria

A análise principal deve priorizar o contraste tumor vs. adjacente e o grau, com ajuste por composição apresentado como análise de atenuação. Estágio, sorologia, Child-Pugh, histologia, fibrose e invasão vascular não sustentam inferência robusta neste desenho. A comparação TCGA × GTEx permanece confundida com estudo/lote e não deve ser interpretada como efeito puro de tumor.
