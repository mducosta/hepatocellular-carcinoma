# 7. Guia de reprodutibilidade

> Arquivo: `docs/reprodutibilidade.md`.

---

## 7.1 Requisitos

- **R ≥ 4.5** (testado com R 4.6.1 no Windows 11).
- Pacotes listados em `dependencias_R.txt`.
- Espaço em disco:
  - ~20 MB para a pasta `analise_covariaveis_deconvolucao/` (matriz processada incluída);
  - ~1,3 GB adicionais se baixar o transcriptoma completo do Toil.
- Memória: ≥ 8 GB recomendado para ler o TPM completo; a matriz filtrada
  (12.401 × 529) usa bem menos.

Verifique a sessão usada:

```bash
cat docs/sessionInfo_analises_exploratorias.txt
```

---

## 7.2 Ordem de execução

```text
etapa_00_baixar_dados_brutos.R        (opcional; baixa dados brutos)
        │
        ▼
etapa_01_preparar_matrizes_e_metadados.R        (metadados + matriz log2 TPM)
        │
        ▼
etapa_02_analise_exploratoria_fenotipos.R  (fenótipos, UMAP, heatmap, marcadores)
        │
        ▼
etapa_03_expressao_diferencial_limma.R (limma em todas as comparações)
        │
        ▼
etapa_04_volcano_plots.R         (12 volcano plots)
        │
        ▼
etapa_05_correcao_lote_ComBat.R   (sensibilidade/limitação do ComBat)
        │
        ▼
etapa_06_fenotipos_clinicos_avancados.R   (missingness, testes, KM, Cox)
        │
        ▼
etapa_07_grupos_alternativos.R   (recidivas, histologia, vírus, TSS)
        │
        ▼
etapa_08_analises_grupos_separados.R      (análises separadas de cada grupo)
        │
        ▼
etapa_10_auditoria_dados_volcano.R  (auditorias)
etapa_11_benchmarks.R                 (benchmarks)
etapa_12_compliance.R                 (checklist final)
        │
        ▼
etapa_09_gerar_relatorio_consolidado.py        (opcional; Markdown consolidado)
```

`executar_todas_as_etapas.R` executa essa sequência automaticamente.

---

## 7.3 Passo a passo

### 7.3.1 Entrar na pasta

```bash
cd "caminho/para/hepatocellular-carcinoma/analise_covariaveis_deconvolucao"
```

### 7.3.2 (Opcional) Baixar dados

Os arquivos pequenos já estão em `dados/raw/`. Para rebaixar:

```bash
Rscript scripts/etapa_00_baixar_dados_brutos.R
```

Para incluir o transcriptoma completo (~1,3 GB):

```bash
BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE Rscript scripts/etapa_00_baixar_dados_brutos.R
```

O script baixa:

- `LIHC_clinicalMatrix` (TCGA, Xena);
- `TCGA-LIHC.clinical.tsv.gz` (GDC, Xena);
- `GTEX_phenotype.gz` e `TcgaTargetGTEX_phenotype.txt.gz` (Toil/Xena);
- `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` (GTEx Portal);
- opcionalmente `TcgaTargetGtex_rsem_gene_tpm.gz`.

### 7.3.3 Preparar dados

Se você baixou o TPM completo:

```bash
TCGA_GTEX_TPM="C:/caminho/TcgaTargetGtex_rsem_gene_tpm.gz" \
  Rscript scripts/etapa_01_preparar_matrizes_e_metadados.R
```

Se o TPM estiver em `$TEMP/hepato_full/TcgaTargetGtex_rsem_gene_tpm.gz`
(caminho usado no desenvolvimento), o script encontra automaticamente.

Se o TPM completo não estiver disponível, o script usa o painel
`../dados/raw/liver.tsv` (213 genes da via LA) e gera
`expressao_LA_pathway.tsv.gz`. As etapas seguintes usam a matriz completa se
ela já existir em `dados/processed/expressao_log2tpm_filtrada.tsv.gz`; caso
contrário, usam o painel LA.

**Saídas esperadas:**

- `dados/processed/metadata_amostras.tsv`
- `dados/processed/metadata_analise.tsv`
- `dados/processed/expressao_log2tpm_filtrada.tsv.gz` (se TPM completo)
- `dados/processed/expressao_LA_pathway.tsv.gz`
- `dados/processed/resumo_preparacao.txt`

### 7.3.4 Análise exploratória

```bash
Rscript scripts/etapa_02_analise_exploratoria_fenotipos.R
```

**Saídas:** `resultados/tabelas/01–12_*.csv`,
`resultados/figuras/01–12_*.png` e `resultados/controle_qualidade/resumo_analise_exploratoria.txt`.

### 7.3.5 Expressão diferencial

```bash
Rscript scripts/etapa_03_expressao_diferencial_limma.R
```

**Saídas:** `resultados/tabelas/expressao_diferencial/expressao_diferencial_*_completo.csv`,
`expressao_diferencial_*_significativos.csv`, `top20_*`, `resumo_expressao_diferencial_todas_comparacoes.csv`,
`resultados/controle_qualidade/resumo_expressao_diferencial.txt`.

> Observação: nessa etapa, o ComBat com covariável biológica falha
> propositalmente (confundimento TCGA × GTEx). A mensagem é capturada e
> registrada; o pipeline continua. A demonstração completa está no script 04.

### 7.3.6 Volcano plots

```bash
Rscript scripts/etapa_04_volcano_plots.R
```

**Saídas:** `resultados/figuras/volcano_01...12...png` (somente PNG) e
`resultados/tabelas/rotulos_volcano_*.csv`.

### 7.3.7 ComBat

```bash
Rscript scripts/etapa_05_correcao_lote_ComBat.R
```

**Saídas:**

- `resultados/tabelas/expressao_diferencial/expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_ComBat_completo.csv`
  (resultado após ComBat sem covariável);
- `resultados/tabelas/expressao_diferencial/comparacao_ComBat_antes_depois.csv`;
- `resultados/figuras/figura_13_UMAP_ComBat_antes_depois.png`;
- `resultados/figuras/figura_14_AFP_ComBat_antes_depois.png`;
- `resultados/controle_qualidade/resumo_ComBat_demonstracao.txt`.

> Esse script sobrescreve `expressao_diferencial_TCGA_Tumor_vs_GTEx_Normal_completo.csv` com a
> mesma análise limma (antes do ComBat), garantindo consistência com o script
> 02.

### 7.3.8 Fenótipos avançados

```bash
Rscript scripts/etapa_06_fenotipos_clinicos_avancados.R
```

**Saídas:** `resultados/tabelas/13–20_*.csv`,
`resultados/figuras/15–25_*.png` e
`resultados/controle_qualidade/resumo_fenotipos_avancados.txt`.

### 7.3.9 Grupos alternativos

```bash
Rscript scripts/etapa_07_grupos_alternativos.R
```

**Saídas:** `resultados/tabelas/21–22_*.csv`,
`resultados/figuras/26–27_*.png` e
`resultados/controle_qualidade/resumo_grupos_alternativos.txt`.

### 7.3.10 Grupos separados

```bash
Rscript scripts/etapa_08_analises_grupos_separados.R
```

**Saídas:** `resultados/tabelas/23–35_*.csv`,
`resultados/figuras/28–35_*.png` e
`resultados/controle_qualidade/resumo_grupos_separados.txt`.

### 7.3.11 Auditorias, benchmarks e compliance

```bash
Rscript scripts/etapa_10_auditoria_dados_volcano.R
Rscript scripts/etapa_11_benchmarks.R
Rscript scripts/etapa_12_compliance.R
```

**Saídas:** `auditorias/auditoria_dados.md`,
`auditorias/auditoria_volcano_ComBat.md`, `auditorias/benchmarks_pipeline.md`,
`auditorias/compliance_pipeline.md` e respectivos CSVs.

### 7.3.12 Relatório Markdown consolidado (opcional)

```bash
python scripts/etapa_09_gerar_relatorio_consolidado.py
```

---

## 7.4 Execução completa

```bash
Rscript scripts/executar_todas_as_etapas.R
```

O `executar_todas_as_etapas.R`:

1. detecta o TPM completo em `TCGA_GTEX_TPM` ou em
   `$TEMP/hepato_full/TcgaTargetGtex_rsem_gene_tpm.gz`;
2. roda 00 → 01 → 02 → 03;
3. roda 04 somente se a matriz completa existir;
4. roda 05 (fenótipos avançados);
5. roda 06 (grupos alternativos);
6. roda 07 (análises separadas dos grupos alternativos);
7. roda 10, 11 e 12 (auditorias, benchmarks e compliance);
8. para em caso de erro.

### 7.4.1 (Opcional) Gerar Markdown consolidado

```bash
python scripts/etapa_09_gerar_relatorio_consolidado.py
```

Não requer dependências externas. Gera `relatorio_consolidado_seminario.md` na raiz,
concatenando os documentos de `docs/` e ajustando os caminhos das figuras.
O relatório é **somente Markdown**; não há versão HTML.

---

## 7.5 Estrutura de saídas

```text
resultados/
├── figuras/
│   ├── figura_01_amostras_por_grupo.png
│   ├── figura_02_sexo_por_grupo.png
│   ├── figura_03_idade_por_grupo.png
│   ├── figura_04_estagio_TCGA.png
│   ├── figura_05_grau_TCGA.png
│   ├── figura_06_estagio_x_grau_TCGA.png
│   ├── figura_07_UMAP_grupo.png
│   ├── figura_08_UMAP_estudo.png
│   ├── figura_09_UMAP_sexo_TCGA.png
│   ├── figura_10_UMAP_estagio_grau_TCGA.png
│   ├── figura_11_heatmap_top50_variaveis.png
│   ├── figura_12_violino_marcadores.png
│   ├── figura_13_UMAP_ComBat_antes_depois.png
│   ├── figura_14_AFP_ComBat_antes_depois.png
│   ├── figura_15_AFP_por_estagio_grau.png
│   ├── figura_16_missingness_clinico.png
│   ├── 17–23_KM_OS_*.png
│   ├── figura_24_forest_cox_multivariado.png
│   ├── figura_25_correlacao_clinica.png
│   ├── figura_26_grupos_dados.png
│   ├── figura_27_histologia_marcadores.png
│   ├── figura_28_recidivas_expressao.png
│   ├── figura_29_recidivas_correlacao.png
│   ├── figura_30_histologia_clinico.png
│   ├── figura_31_viral_clinico.png
│   ├── figura_32_viral_marcadores.png
│   ├── figura_33_tss_distribuicao.png
│   ├── figura_34_tss_umap.png
│   ├── figura_35_tss_marcadores_heatmap.png
│   └── volcano_01...12...png
├── tabelas/
│   ├── 01–35_*.csv
│   └── DE/
│       ├── expressao_diferencial_*_completo.csv
│       ├── expressao_diferencial_*_significativos.csv
│       ├── top20_*.csv
│       ├── resumo_expressao_diferencial_todas_comparacoes.csv
│       └── comparacao_ComBat_antes_depois.csv
└── controle_qualidade/
    ├── log_etapa_01_preparar_matrizes_e_metadados.txt
    ├── log_etapa_02_analise_exploratoria.txt
    ├── log_etapa_03_expressao_diferencial.txt
    ├── log_etapa_04_volcano.txt
    ├── log_etapa_05_ComBat.txt
    ├── log_etapa_06_fenotipos_avancados.txt
    ├── log_etapa_07_grupos_alternativos.txt
    ├── resumo_grupos_separados.txt
    ├── log_etapa_10_auditoria.txt
    ├── log_etapa_11_benchmarks.txt
    ├── log_etapa_12_compliance.txt
    ├── resumo_analise_exploratoria.txt
    ├── resumo_expressao_diferencial.txt
    ├── resumo_ComBat_demonstracao.txt
    ├── coordenadas_UMAP_ComBat.csv
    └── resumo_fenotipos_avancados.txt
```

Os relatórios de auditoria, benchmark e compliance ficam em `auditorias/`:

```text
auditorias/
├── auditoria_dados.md
├── auditoria_volcano_ComBat.md
├── benchmarks_pipeline.md
├── compliance_pipeline.md
├── auditoria_dados_resumo.csv
├── auditoria_volcano_ComBat_diagnostico.csv
├── benchmarks_pipeline_tabela.csv
└── compliance_pipeline_tabela.csv
```

---

## 7.6 Solução de problemas

### “O arquivo TPM não é encontrado”

- Verifique se `TCGA_GTEX_TPM` aponta para o arquivo `.gz`.
- No Windows, use barras normais (`C:/...`) ou barras invertidas escapadas.
- Se não tiver o TPM completo, rode sem ele: o script 00 usa o painel LA. Se a
  matriz completa já existir em `dados/processed/expressao_log2tpm_filtrada.tsv.gz`,
  as etapas seguintes continuam usando essa matriz.

### “The covariate is confounded with batch!”

- É o comportamento **esperado** no script 04 para `mod = ~ grupo`.
- O script captura o erro e segue para a demonstração sem covariável.
- Não é um bug; é a limitação central do ComBat neste desenho.

### “Poucas amostras” ou “NA”

- Confira `metadata_analise.tsv`; amostras sem idade/sexo podem ser removidas
  por `complete.cases` nas comparações ajustadas.
- O TCGA tem 1 tumor sem idade; o GTEx tem idade para todos os 110.

### Caminhos e acentuação

- A pasta `analise_covariaveis_deconvolucao` usa apenas caracteres ASCII, sem
  acentos ou espaços, e funciona em qualquer locale. Os scripts usam caminhos
  relativos e `normalizePath`, então funcionam no Windows 11 com locale pt-BR.

### Diferenças numéricas

- Pequenas diferenças podem ocorrer entre versões de R/Bioconductor.
- As versões usadas estão em `docs/sessionInfo_analises_exploratorias.txt`.
- O número de genes após o filtro pode variar se o mapeamento
  Ensembl → símbolo mudar.

---

## 7.7 Integridade e versionamento

- A pasta contém apenas resultados leves (figuras, tabelas, matriz processada).
- O TPM bruto de 1,3 GB **não** deve ser versionado.
- Se for usar Git, adicione ao `.gitignore`:

```gitignore
analise_covariaveis_deconvolucao/dados/raw/TcgaTargetGtex_rsem_gene_tpm.gz
analise_covariaveis_deconvolucao/resultados/controle_qualidade/*.log
```

- Os scripts são determinísticos, exceto:
  - `ComBat` (usa `set.seed(2026)` no script 04);
  - `ggrepel` (posicionamento de rótulos; `seed = 2026` no volcano);
  - pequenas variações de versão de pacote.
