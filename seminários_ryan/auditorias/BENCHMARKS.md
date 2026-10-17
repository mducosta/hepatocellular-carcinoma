# Benchmarks computacionais e de concordância

**Data:** 2026-10-17

## 1. Tempo de execução

|   |benchmark                             | valor|unidade |observacao                                    |
|:--|:-------------------------------------|-----:|:-------|:---------------------------------------------|
|1  |Carga de metadados + matriz           |  1.50|s       |12401 genes x 529 amostras                    |
|3  |UMAP (2.000 genes, 529 amostras)      |  4.34|s       |n_neighbors=30, min_dist=0.25, cosine, pca=50 |
|4  |limma com covariáveis (tumor vs GTEx) |  0.69|s       |12401 genes x 478 amostras                    |
|5  |limma sem covariáveis (tumor vs GTEx) |  0.26|s       |mesmo contraste, sem sexo/idade               |
|6  |limma com trend = FALSE               |  0.47|s       |sensibilidade do modelo de variância          |

## 2. Concordância estatística

|   |benchmark                                    |     valor|unidade |observacao                     |
|:--|:--------------------------------------------|---------:|:-------|:------------------------------|
|2  |Memória da matriz                            |    53.300|MB      |matriz numérica em memória     |
|7  |DEGs com covariáveis                         |  3388.000|genes   |FDR<0.05 e &#124;logFC&#124;>1 |
|8  |DEGs sem covariáveis                         |  3316.000|genes   |FDR<0.05 e &#124;logFC&#124;>1 |
|9  |Jaccard com vs sem covariáveis               |     0.948|índice  |interseção/união dos DEGs      |
|10 |Correlação de logFC (com vs sem covariáveis) |     0.999|rho     |Spearman                       |
|11 |DEGs trend = TRUE                            |  3388.000|genes   |                               |
|12 |DEGs trend = FALSE                           |  3389.000|genes   |                               |
|13 |Jaccard trend TRUE vs FALSE                  |     1.000|índice  |                               |
|14 |Genes significativos apenas por FDR          | 10109.000|genes   |FDR<0.05                       |
|15 |Genes significativos por FDR + logFC         |  3388.000|genes   |FDR<0.05 e &#124;logFC&#124;>1 |
|16 |Jaccard FDR vs FDR+logFC                     |     0.335|índice  |                               |
|17 |DEGs tumor vs adjacente (não pareado)        |  4130.000|genes   |                               |
|18 |DEGs tumor vs adjacente (pareado)            |  3710.000|genes   |                               |
|19 |Jaccard pareado vs não pareado               |     0.764|índice  |                               |
|20 |Correlação de logFC (pareado vs não pareado) |     0.969|rho     |Spearman                       |
|21 |DEGs antes do ComBat                         |  3388.000|genes   |                               |
|22 |DEGs depois do ComBat (sem covariável)       |     1.000|genes   |                               |
|23 |Jaccard ComBat antes vs depois               |     0.000|índice  |                               |
|24 |Correlação de logFC (ComBat antes vs depois) |     0.647|rho     |Spearman                       |

## 3. Interpretação

- O **limma** é rápido o suficiente para 12.401 genes × 478 amostras; o
  custo dominante é a leitura da matriz e a geração das figuras.
- O **UMAP** é viável em poucos segundos com 2.000 genes e 529 amostras,
  usando `n_threads = 1` para maior determinismo.
- A inclusão de **sexo e idade** altera parte dos DEGs; a correlação de
  logFC entre modelos costuma ser alta, mas o conjunto de genes
  significativos não é idêntico. Isso justifica reportar o modelo ajustado
  como principal.
- A comparação **pareado vs não pareado** mostra alta concordância, o que
  reforça a robustez do sinal tumor vs. adjacente.
- A demonstração do **ComBat** reduz drasticamente os DEGs; o Jaccard baixo
  e a correlação de logFC baixa confirmam a sobrecorreção.
- O critério **FDR + |logFC| > 1** é mais conservador que FDR apenas; os
  genes excluídos pelo corte de logFC podem ser consultados nas tabelas
  completas.

## 4. Arquivos

- `auditorias/BENCHMARKS_tabela.csv`
- `resultados/tabelas/DE/`

## 5. Sessão R

```text
R version 4.6.1 (2026-06-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26340)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Portuguese_Brazil.utf8  LC_CTYPE=Portuguese_Brazil.utf8   
[3] LC_MONETARY=Portuguese_Brazil.utf8 LC_NUMERIC=C                      
[5] LC_TIME=Portuguese_Brazil.utf8    

time zone: America/Sao_Paulo
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] knitr_1.52        limma_3.68.0      data.table_1.18.4

loaded via a namespace (and not attached):
 [1] RcppAnnoy_0.0.23  codetools_0.2-20  xfun_0.62         Matrix_1.7-5     
 [5] lattice_0.22-9    splines_4.6.1     R.utils_2.13.0    cli_3.6.6        
 [9] R.methodsS3_1.8.2 grid_4.6.1        statmod_1.5.1     compiler_4.6.1   
[13] R.oo_1.27.1       tools_4.6.1       RSpectra_0.16-2   irlba_2.4.1      
[17] uwot_0.2.4        evaluate_1.0.5    Rcpp_1.1.2        otel_0.2.0       
[21] rlang_1.3.0      
```
