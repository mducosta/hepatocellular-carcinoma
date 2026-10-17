# Compliance e integridade dos artefatos

**Data:** 2026-10-17

## 1. Checklist

|item                                | status |evidencia                                                                                                                                                                                                                                                                                                    |
|:-----------------------------------|:------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|Arquivos obrigatórios               |   OK   |ausentes: 0                                                                                                                                                                                                                                                                                                  |
|Figuras sem PDF                     |   OK   |PDFs encontrados: 0                                                                                                                                                                                                                                                                                          |
|PNGs válidos e não vazios           |   OK   |PNGs: 37; suspeitos de vazio: 0                                                                                                                                                                                                                                                                              |
|Volcano_12 não vazia                |   OK   |fração não branca: 0.0948                                                                                                                                                                                                                                                                                    |
|Tabelas CSV válidas                 |   OK   |CSVs: 73; vazios/ilegíveis: 0                                                                                                                                                                                                                                                                                |
|Colunas críticas sem NA             |   OK   |resultados/tabelas/DE/DEG_TCGA_Tumor_vs_GTEx_Normal_completo.csv: gene=0; logFC=0; P.Value=0; adj.P.Val=0 &#124; resultados/tabelas/14_fenotipos_testes.csv: variavel_1=0; variavel_2=0; p_valor=0 &#124; resultados/tabelas/18_cox_multivariado.csv: term=0; estimate=0; conf.low=0; conf.high=0; p.value=0 |
|Scripts R parseáveis                |   OK   |scripts: 11; erros: 0                                                                                                                                                                                                                                                                                        |
|Sementes de reprodutibilidade       |   OK   |01_analise_exploratoria.R=TRUE; 03_volcano_plots.R=TRUE; 04_ComBat_demonstracao.R=TRUE; 05_fenotipos_avancados.R=TRUE                                                                                                                                                                                        |
|Imagens dos docs existem            |   OK   |docs: 9; imagens ausentes: 0                                                                                                                                                                                                                                                                                 |
|Imagens do HTML existem             |   OK   |ausentes: 0                                                                                                                                                                                                                                                                                                  |
|Docs sem PDF                        |   OK   |nenhuma menção                                                                                                                                                                                                                                                                                               |
|Docs sem PCA (substituído por UMAP) |   OK   |nenhuma menção                                                                                                                                                                                                                                                                                               |
|Colunas do metadado                 |   OK   |ausentes: 0                                                                                                                                                                                                                                                                                                  |
|Contagem de amostras                |   OK   |n = 529 (esperado 529)                                                                                                                                                                                                                                                                                       |
|Sem recidivas na análise            |   OK   |amostras -02: 0                                                                                                                                                                                                                                                                                              |

## 2. Critérios avaliados

- existência dos scripts, tabelas, figuras e documentos essenciais;
- ausência de figuras em PDF (apenas PNG);
- PNGs válidos, legíveis e não vazios;
- `Volcano_12` com conteúdo visual mínimo;
- CSVs legíveis, não vazios e com colunas críticas sem NA;
- scripts R parseáveis e com sementes de reprodutibilidade;
- imagens referenciadas nos docs e no HTML existentes;
- ausência de menções a PDF e a PCA nos documentos;
- metadados com colunas esperadas, 529 amostras e sem recidivas;
- versão do R e pacotes registrada em `docs/sessionInfo.txt`.

## 3. Resultado

**Todos os itens passaram.**

## 4. Observações

- A checagem de “não vazio” usa a fração de pixels não brancos; ela detecta
  figuras totalmente brancas, mas não substitui inspeção visual.
- A checagem de PCA/PDF é textual e complementa a inspeção dos scripts.
- A ausência de validação experimental e o confundimento TCGA × GTEx são
  limitações científicas, não falhas de compliance.

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
[1] png_0.1-9         knitr_1.52        data.table_1.18.4

loaded via a namespace (and not attached):
[1] compiler_4.6.1 cli_3.6.6      tools_4.6.1    otel_0.2.0     xfun_0.62     
[6] rlang_1.3.0    evaluate_1.0.5
```
