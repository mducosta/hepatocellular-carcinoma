# Auditoria da figura `volcano_12_TCGA_vs_GTEx_ComBat.png`

**Data da auditoria:** 2026-10-17

## 1. Queixa

A figura `volcano_12` foi reportada como “sem nada” (visualmente vazia).

## 2. Diagnóstico

|metrica                                |valor           |
|:--------------------------------------|:---------------|
|Arquivo existe                         |TRUE            |
|Dimensões (px)                         |1800 x 3600 x 3 |
|Fração de pixels não brancos           |0.0948          |
|Fração de pixels avermelhados          |0.00980         |
|Fração de pixels azulados              |0.00804         |
|Genes testados (depois do ComBat)      |12401           |
|DEGs significativos (depois do ComBat) |1               |
|DEGs significativos (antes do ComBat)  |3388            |
|Menor FDR (depois)                     |7.37e-21        |
|Maior &#124;logFC&#124; (depois)       |6.149           |
|Gene significativo (depois)            |LINC01419       |

**A figura não estava tecnicamente em branco.** A auditoria de pixels
mostra que a imagem original tinha aproximadamente **3,65%** de pixels não
brancos (eixos, título, legenda e poucos pontos), contra **7,39%** do
`volcano_03` (comparação tumor vs. normal antes do ComBat). O problema era
**científico/visual**: após o ComBat sem covariável biológica, restou apenas
**1 gene significativo** (LINC01419). Os outros 12.400 genes ficaram
amontoados perto da origem (logFC ≈ 0 e FDR ≈ 1), de modo que o painel
parecia vazio.

## 3. Causa raiz

1. O ComBat foi aplicado com `batch = estudo` (TCGA/GTEx) e **sem**
   covariável de condição (`mod = NULL`).
2. Como estudo e condição são perfeitamente confundidos, a correção removeu
   também a diferença tumor-normal.
3. O resultado esperado é a quase ausência de DEGs — **não** um erro de
   renderização do `ggplot2`.
4. A tentativa “correta” (`mod = ~ grupo`) é recusada pelo próprio ComBat:
   `The covariate is confounded with batch!`.

## 4. Correção aplicada

- A figura foi refeita como **dois painéis**: antes do ComBat (3.388 DEGs)
  e depois do ComBat (1 DEG).
- O painel “depois” agora tem escala livre e rótulo do único gene
  significativo (LINC01419), além de subtítulo com a contagem de DEGs.
- A figura corrigida tem **9,48%** de pixels não brancos, quase 3× a
  original, e comunica explicitamente a sobrecorreção.

## 5. Evidências

- Dados: `resultados/tabelas/expressao_diferencial/expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_ComBat_completo.csv`
- Comparação: `resultados/tabelas/expressao_diferencial/comparacao_ComBat_antes_depois.csv`
- Figura corrigida: `resultados/figuras/volcano_12_TCGA_vs_GTEx_ComBat.png`
- Log: `resultados/controle_qualidade/log_etapa_05_ComBat.txt`

## 6. Conclusão

A figura original era um **artefato científico esperado** da sobrecorreção
do ComBat, não um bug de código. A versão corrigida é auditável, mostra o
antes/depois e deixa explícito que o resultado “vazio” é a demonstração da
limitação, não um resultado biológico.
