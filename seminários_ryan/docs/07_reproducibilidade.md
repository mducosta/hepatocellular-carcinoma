# 7. Guia de reprodutibilidade

> Arquivo: `docs/07_reproducibilidade.md`.

---

## 7.1 Requisitos

- **R ≥ 4.5** (testado com R 4.6.1 no Windows 11).
- Pacotes listados em `requirements_R.txt`.
- Espaço em disco:
  - ~20 MB para a pasta `seminários_ryan/` (matriz processada incluída);
  - ~1,3 GB adicionais se baixar o transcriptoma completo do Toil.
- Memória: ≥ 8 GB recomendado para ler o TPM completo; a matriz filtrada
  (12.401 × 529) usa bem menos.

Verifique a sessão usada:

```bash
cat docs/sessionInfo.txt
```

---

## 7.2 Ordem de execução

```text
00_download_dados.R        (opcional; baixa dados brutos)
        │
        ▼
00_preparar_dados.R        (metadados + matriz log2 TPM)
        │
        ▼
01_analise_exploratoria.R  (fenótipos, UMAP, heatmap, marcadores)
        │
        ▼
02_expressao_diferencial.R (limma em todas as comparações)
        │
        ▼
03_volcano_plots.R         (12 volcano plots)
        │
        ▼
04_ComBat_demonstracao.R   (sensibilidade/limitação do ComBat)
        │
        ▼
05_fenotipos_avancados.R   (missingness, testes, KM, Cox)
        │
        ▼
10_auditoria_dados_e_volcano.R  (auditorias)
11_benchmarks.R                 (benchmarks)
12_compliance.R                 (checklist final)
        │
        ▼
09_gerar_relatorio_html.py      (opcional; HTML consolidado)
```

`run_all.R` executa essa sequência automaticamente.

---

## 7.3 Passo a passo

### 7.3.1 Entrar na pasta

```bash
cd "caminho/para/hepatocellular-carcinoma/seminários_ryan"
```

### 7.3.2 (Opcional) Baixar dados

Os arquivos pequenos já estão em `dados/raw/`. Para rebaixar:

```bash
Rscript scripts/00_download_dados.R
```

Para incluir o transcriptoma completo (~1,3 GB):

```bash
BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE Rscript scripts/00_download_dados.R
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
  Rscript scripts/00_preparar_dados.R
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
Rscript scripts/01_analise_exploratoria.R
```

**Saídas:** `resultados/tabelas/01–12_*.csv`,
`resultados/figuras/01–12_*.png` e `resultados/qc/resumo_analise_exploratoria.txt`.

### 7.3.5 Expressão diferencial

```bash
Rscript scripts/02_expressao_diferencial.R
```

**Saídas:** `resultados/tabelas/DE/DEG_*_completo.csv`,
`DEG_*_significativos.csv`, `TOP20_*`, `RESUMO_DEGs_TODAS_COMPARACOES.csv`,
`resultados/qc/resumo_expressao_diferencial.txt`.

> Observação: nessa etapa, o ComBat com covariável biológica falha
> propositalmente (confundimento TCGA × GTEx). A mensagem é capturada e
> registrada; o pipeline continua. A demonstração completa está no script 04.

### 7.3.6 Volcano plots

```bash
Rscript scripts/03_volcano_plots.R
```

**Saídas:** `resultados/figuras/Volcano_01...12...png` (somente PNG) e
`resultados/tabelas/ROTULOS_Volcano_*.csv`.

### 7.3.7 ComBat

```bash
Rscript scripts/04_ComBat_demonstracao.R
```

**Saídas:**

- `resultados/tabelas/DE/DEG_TCGA_Tumor_vs_GTEx_Normal_ComBat_completo.csv`
  (resultado após ComBat sem covariável);
- `resultados/tabelas/DE/COMPARACAO_ComBat_antes_depois.csv`;
- `resultados/figuras/13_UMAP_ComBat_antes_depois.png`;
- `resultados/figuras/14_AFP_ComBat_antes_depois.png`;
- `resultados/qc/resumo_ComBat_demonstracao.txt`.

> Esse script sobrescreve `DEG_TCGA_Tumor_vs_GTEx_Normal_completo.csv` com a
> mesma análise limma (antes do ComBat), garantindo consistência com o script
> 02.

### 7.3.8 Fenótipos avançados

```bash
Rscript scripts/05_fenotipos_avancados.R
```

**Saídas:** `resultados/tabelas/13–20_*.csv`,
`resultados/figuras/15–25_*.png` e
`resultados/qc/resumo_fenotipos_avancados.txt`.

### 7.3.9 Auditorias, benchmarks e compliance

```bash
Rscript scripts/10_auditoria_dados_e_volcano.R
Rscript scripts/11_benchmarks.R
Rscript scripts/12_compliance.R
```

**Saídas:** `auditorias/AUDITORIA_DADOS.md`,
`auditorias/AUDITORIA_VOLCANO_12.md`, `auditorias/BENCHMARKS.md`,
`auditorias/COMPLIANCE.md` e respectivos CSVs.

### 7.3.10 Relatório HTML consolidado (opcional)

```bash
python scripts/09_gerar_relatorio_html.py
```

---

## 7.4 Execução completa

```bash
Rscript scripts/run_all.R
```

O `run_all.R`:

1. detecta o TPM completo em `TCGA_GTEX_TPM` ou em
   `$TEMP/hepato_full/TcgaTargetGtex_rsem_gene_tpm.gz`;
2. roda 00 → 01 → 02 → 03;
3. roda 04 somente se a matriz completa existir;
4. roda 05 (fenótipos avançados);
5. roda 10, 11 e 12 (auditorias, benchmarks e compliance);
6. para em caso de erro.

### 7.4.1 (Opcional) Gerar HTML consolidado

```bash
python scripts/09_gerar_relatorio_html.py
```

Requer `pip install markdown`. Gera `relatorio_completo.html` na raiz,
concatenando os documentos de `docs/` e as figuras.

---

## 7.5 Estrutura de saídas

```text
resultados/
├── figuras/
│   ├── 01_amostras_por_grupo.png
│   ├── 02_sexo_por_grupo.png
│   ├── 03_idade_por_grupo.png
│   ├── 04_estagio_TCGA.png
│   ├── 05_grau_TCGA.png
│   ├── 06_estagio_x_grau_TCGA.png
│   ├── 07_UMAP_grupo.png
│   ├── 08_UMAP_estudo.png
│   ├── 09_UMAP_sexo_TCGA.png
│   ├── 10_UMAP_estagio_grau_TCGA.png
│   ├── 11_heatmap_top50_variaveis.png
│   ├── 12_boxplot_marcadores.png
│   ├── 13_UMAP_ComBat_antes_depois.png
│   ├── 14_AFP_ComBat_antes_depois.png
│   ├── 15_AFP_por_estagio_grau.png
│   ├── 16_missingness_clinico.png
│   ├── 17–23_KM_OS_*.png
│   ├── 24_forest_cox_multivariado.png
│   ├── 25_correlacao_clinica.png
│   └── Volcano_01...12...png
├── tabelas/
│   ├── 01–20_*.csv
│   └── DE/
│       ├── DEG_*_completo.csv
│       ├── DEG_*_significativos.csv
│       ├── TOP20_*.csv
│       ├── RESUMO_DEGs_TODAS_COMPARACOES.csv
│       └── COMPARACAO_ComBat_antes_depois.csv
└── qc/
    ├── log_00_preparar_dados.txt
    ├── log_01_analise_exploratoria.txt
    ├── log_02_expressao_diferencial.txt
    ├── log_03_volcano.txt
    ├── log_04_ComBat.txt
    ├── log_05_fenotipos_avancados.txt
    ├── log_10_auditoria.txt
    ├── log_11_benchmarks.txt
    ├── log_12_compliance.txt
    ├── resumo_analise_exploratoria.txt
    ├── resumo_expressao_diferencial.txt
    ├── resumo_ComBat_demonstracao.txt
    ├── UMAP_ComBat_coordenadas.csv
    └── resumo_fenotipos_avancados.txt
```

Os relatórios de auditoria, benchmark e compliance ficam em `auditorias/`:

```text
auditorias/
├── AUDITORIA_DADOS.md
├── AUDITORIA_VOLCANO_12.md
├── BENCHMARKS.md
├── COMPLIANCE.md
├── auditoria_dados_resumo.csv
├── auditoria_volcano12_diagnostico.csv
├── BENCHMARKS_tabela.csv
└── COMPLIANCE_tabela.csv
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

### Acentuação em caminhos

- A pasta `seminários_ryan` tem acento. Em sistemas com locale não UTF-8,
  renomeie para `seminarios_ryan` ou execute a partir da pasta com
  `cd`. Os scripts usam caminhos relativos e `normalizePath`, então funcionam
  no Windows 11 com locale pt-BR.

### Diferenças numéricas

- Pequenas diferenças podem ocorrer entre versões de R/Bioconductor.
- As versões usadas estão em `docs/sessionInfo.txt`.
- O número de genes após o filtro pode variar se o mapeamento
  Ensembl → símbolo mudar.

---

## 7.7 Integridade e versionamento

- A pasta contém apenas resultados leves (figuras, tabelas, matriz processada).
- O TPM bruto de 1,3 GB **não** deve ser versionado.
- Se for usar Git, adicione ao `.gitignore`:

```gitignore
seminários_ryan/dados/raw/TcgaTargetGtex_rsem_gene_tpm.gz
seminários_ryan/resultados/qc/*.log
```

- Os scripts são determinísticos, exceto:
  - `ComBat` (usa `set.seed(2026)` no script 04);
  - `ggrepel` (posicionamento de rótulos; `seed = 2026` no volcano);
  - pequenas variações de versão de pacote.
