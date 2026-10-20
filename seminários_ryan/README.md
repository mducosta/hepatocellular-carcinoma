# Seminário Ryan — TCGA-LIHC × GTEx (fígado)

Análise exploratória, clínica e de **expressão diferencial** do
**hepatocarcinoma (TCGA-LIHC)** em comparação com **fígado normal (GTEx)**,
com explicações didáticas sobre **ComBat**, **logFC no limma** e
**volcano plots**.

> Pasta: `seminários_ryan/` (dentro do repositório `hepatocellular-carcinoma`).
> Todos os scripts, resultados, figuras e explicações estão organizados aqui.

---

## Sumário

1. [Pergunta e desenho](#1-pergunta-e-desenho)
2. [Dados](#2-dados)
3. [Estrutura da pasta](#3-estrutura-da-pasta)
4. [Resultados principais](#4-resultados-principais)
5. [Documentação didática](#5-documentação-didática)
6. [Como reproduzir](#6-como-reproduzir)
7. [Limitações](#7-limitações)
8. [Declaração de uso de IA](#8-declaração-de-uso-de-inteligência-artificial)
9. [Referências](#9-referências)

---

## 1. Pergunta e desenho

**Pergunta:** como o transcriptoma do hepatocarcinoma difere do fígado normal
e de tecido hepático adjacente, e como fatores clínicos (sexo, idade, estágio,
grau) se associam a esses padrões?

**Comparações realizadas:**

| Bloco | Comparação |
|---|---|
| TCGA vs. TCGA | Tumor primário vs. tecido adjacente (não pareado e pareado) |
| TCGA vs. TCGA | Estágio III–IV vs. I–II |
| TCGA vs. TCGA | Grau G3–G4 vs. G1–G2 |
| TCGA vs. TCGA | Masculino vs. feminino; idade ≥65 vs. <65 |
| GTEx vs. GTEx | Masculino vs. feminino |
| GTEx vs. GTEx | Idade ≥60 vs. <60 |
| TCGA vs. GTEx | Tumor primário vs. fígado normal |
| TCGA vs. GTEx | Tecido adjacente vs. fígado normal |
| Via KEGG hsa05417 | Tumor vs. normal e tumor vs. adjacente |
| Sensibilidade | ComBat (demonstração de sobrecorreção) |

---

## 2. Dados

| Grupo | n |
|---|---:|
| TCGA-LIHC — tumor primário | 369 |
| TCGA-LIHC — tecido adjacente | 50 |
| GTEx — fígado normal | 110 |
| **Total (sem recidivas)** | **529** |

- **Expressão:** UCSC Xena / Toil, `TcgaTargetGtex_rsem_gene_tpm`,
  log2(TPM + 0,001).
- **Genes após filtro:** 12.401.
- **Painel KEGG hsa05417:** 213 genes.
- **Clínico:** TCGA GDC/`LIHC_clinicalMatrix`; GTEx v8 (sexo e faixa etária).
- **Excluídos:** 2 tumores recidivados (`-02`).

> O transcriptoma completo bruto (~1,3 GB) não está versionado. A matriz
> processada (17 MB) está em
> `dados/processed/expressao_log2tpm_filtrada.tsv.gz`.

---

## 3. Estrutura da pasta

```text
seminários_ryan/
├── README.md
├── relatorio_completo.md            # versão Markdown consolidada dos docs
├── requirements_R.txt
├── scripts/
│   ├── 00_download_dados.R          # download dos dados brutos
│   ├── 00_preparar_dados.R          # metadados + matriz de expressão
│   ├── 01_analise_exploratoria.R    # fenótipos, UMAP, heatmap, marcadores
│   ├── 02_expressao_diferencial.R   # limma: todas as comparações
│   ├── 03_volcano_plots.R           # volcano plots (PNG)
│   ├── 04_ComBat_demonstracao.R     # demonstração e limitação do ComBat
│   ├── 05_fenotipos_avancados.R     # missingness, testes, KM, Cox
│   ├── 06_grupos_alternativos.R     # inventário: recidivas, histologia, vírus, TSS
│   ├── 07_grupos_separados.R        # análises separadas de cada grupo
│   ├── 10_auditoria_dados_e_volcano.R
│   ├── 11_benchmarks.R
│   ├── 12_compliance.R
│   ├── 09_gerar_relatorio_md.py     # (opcional) Markdown consolidado
│   └── run_all.R                    # executa tudo em ordem
├── dados/
│   ├── raw/                         # clínico/fenótipo (pequenos)
│   └── processed/
│       ├── metadata_amostras.tsv
│       ├── metadata_analise.tsv
│       ├── expressao_log2tpm_filtrada.tsv.gz
│       ├── expressao_LA_pathway.tsv.gz
│       └── resumo_preparacao.txt
├── resultados/
│   ├── tabelas/
│   │   ├── 01–35_*.csv              # fenótipos, UMAP, testes, Cox, grupos
│   │   └── DE/                      # resultados de expressão diferencial
│   ├── figuras/                     # UMAP, heatmap, KM, volcanos (PNG)
│   └── qc/                          # logs e resumos de execução
├── auditorias/
│   ├── AUDITORIA_DADOS.md
│   ├── AUDITORIA_VOLCANO_12.md
│   ├── BENCHMARKS.md
│   └── COMPLIANCE.md
└── docs/
    ├── 01_analise_exploratoria.md
    ├── 02_ComBat.md
    ├── 03_logFC_limma.md
    ├── 04_volcano_plots.md
    ├── 05_relatorio_final.md
    ├── 06_dicionario_dados.md
    ├── 07_reproducibilidade.md
    ├── 08_fenotipos_avancados.md
    ├── 09_auditorias_benchmarks.md
    ├── 10_grupos_alternativos.md
    ├── 11_grupos_separados.md
    └── sessionInfo.txt
```

---

## 4. Resultados principais

### 4.1 Fenótipos

- **Sexo:** 120 mulheres / 249 homens no TCGA tumor; 33/77 no GTEx normal;
  sem diferença significativa (p = 0,703).
- **Idade:** TCGA tumor 59,4 ± 13,5 anos vs. GTEx 53,7 ± 11,1 anos
  (Wilcoxon p = 3,96 × 10⁻⁶).
- **Estágio (TCGA tumor):** I = 169, II = 86, III = 85, IV = 5.
- **Grau (TCGA tumor):** G1 = 55, G2 = 177, G3 = 121, G4 = 11.
- **UMAP:** UMAP2 associa-se a grau (ρ = −0,346), estágio (ρ = −0,222)
  e idade (ρ = 0,207); o UMAP não fornece percentual de variância explicada.

### 4.2 Expressão diferencial (FDR < 0,05 e |logFC| > 1)

| Comparação | Up | Down | Total |
|---|---:|---:|---:|
| TCGA tumor vs. adjacente | 3.322 | 808 | 4.130 |
| TCGA tumor vs. adjacente (pareado) | 2.934 | 776 | 3.710 |
| TCGA tumor vs. GTEx normal | 2.245 | 1.143 | 3.388 |
| TCGA adjacente vs. GTEx normal | 584 | 1.966 | 2.550 |
| TCGA masculino vs. feminino | 155 | 64 | 219 |
| TCGA G3–G4 vs. G1–G2 | 61 | 111 | 172 |
| GTEx ≥60 vs. <60 | 2 | 113 | 115 |
| TCGA ≥65 vs. <65 | 43 | 30 | 73 |
| TCGA estágio III–IV vs. I–II | 3 | 43 | 46 |
| Via LA: TCGA tumor vs. GTEx | 14 | 32 | 46 |
| Via LA: TCGA tumor vs. adjacente | 11 | 33 | 44 |
| GTEx masculino vs. feminino | 8 | 5 | 13 |

### 4.3 Fenótipos avançados

- **Missingness (TCGA tumor):** invasão vascular 70,7%, sorologia viral 55,8%,
  fibrose 43,1%, Child-Pugh 35,8% e AFP 25,2%.
- **Associações significativas (BH < 0,05):** Child-Pugh × estágio
  (p = 5,5 × 10⁻⁸), AFP × grau (p = 4,3 × 10⁻⁶), vírus × grau/estágio
  (p ≈ 1,1–1,4 × 10⁻²), fibrose × estágio (p = 6,7 × 10⁻³), entre outras.
- **Sobrevida:** apenas estágio foi significativo no log-rank
  (p = 5,18 × 10⁻⁷); grau, sexo, idade, AFP e fibrose não.
- **Cox multivariado:** estágio III–IV vs. I–II HR = 2,61 (1,79–3,81;
  p = 6,0 × 10⁻⁷). No Modelo B, idade HR = 1,03/ano (p = 0,0165); AFP e
  fibrose não foram independentes.
- **Riscos proporcionais:** `cox.zph` global p = 0,386 (Modelo A) e 0,153
  (Modelo B).

### 4.4 Outros grupos nos dados

- **Recidivas:** 2 amostras `Recurrent Tumor`, excluídas da análise principal.
- **Histologia:** HCC clássico 359; hepatocolangiocarcinoma misto 7;
  fibrolamelar 3.
- **Vírus:** HBV+HCV 85; HBV 59; HCV 19; unknown 206.
- **TSS:** 38 centros TCGA; `DD` concentra 144 tumores (39%).
- **Marcadores por histologia:** o grupo misto tem KRT19 e EPCAM muito altos,
  consistente com componente colangiocelular; fibrolamelar tem AFP/GPC3 mais
  baixos. n pequeno: apenas descritivo.
- **Análises separadas:** clínica das 2 recidivas + expressão no painel LA;
  clínica e marcadores por subtipo; clínica, testes e marcadores por grupo
  viral; distribuição, clínica, marcadores e UMAP por TSS.
- **Vírus:** `Unknown` difere de casos conhecidos em estágio, grau e fibrose;
  a missingness é informativa e não deve ser tratada como “sem infecção”.
- **TSS:** centros com n ≥ 10 diferem em clínica e marcadores; o TSS explica
  ~8% de UMAP1 e ~7% de UMAP2, compatível com efeito de centro/lote.
- Detalhes em `docs/10_grupos_alternativos.md` e
  `docs/11_grupos_separados.md`.

### 4.5 Achados biológicos

- **Perda de identidade hepática:** CLEC4M, CLEC4G, CLEC1B, FCN2/3, CYP1A2,
  GDF2, MARCO, CRHBP, ADAMTS13.
- **Proliferação:** BIRC5, UBE2T, CDC25C, AURKB, CCNB1, CDC20, NUF2, KIF2C.
- **Remodelamento vascular/matricial:** PLVAP, CD34, MMP9, MMP1, VIM, SPP1,
  COL1A1, VCAM1.
- **Hipóxia/estágio avançado:** CA9, MIR210HG, SLC6A8 (up); GNMT, CYP2C9,
  ACOT12, FTCD, HPX, SLC10A1 (down).
- **Controle positivo:** genes do cromossomo Y e XIST na comparação por sexo.

### 4.6 ComBat

A comparação TCGA × GTEx é **estruturalmente confundida**: todo tumor é TCGA e
todo normal é GTEx. O ComBat com `mod = ~ grupo` recusa o modelo
(`The covariate is confounded with batch!`). Sem covariável, ele reduz
**3.388 DEGs para 1**, demonstrando **sobrecorreção**. Por isso, o ComBat não
foi aplicado à análise principal. Detalhes em `docs/02_ComBat.md`.

### 4.7 Auditorias, benchmarks e compliance

- **Auditoria de dados:** 22 verificações; único alerta = confundimento
  estrutural TCGA × GTEx (Cramér's V = 1).
- **Auditoria do Volcano_12:** a figura original parecia vazia porque restou
  apenas 1 DEG após o ComBat; foi refeita como painel antes/depois
  (3.388 → 1 DEG) e agora tem ≈9,48% de pixels não brancos.
- **Benchmarks:** UMAP ~4,3 s; limma ~0,7 s (12.401 genes × 478 amostras);
  Jaccard pareado vs. não pareado = 0,764; ComBat antes vs. depois = 0.
- **Compliance:** todos os itens do checklist passaram (sem PDFs, PNGs
  válidos, scripts parseáveis, sementes, links e metadados).

---

## 5. Documentação didática

| Documento | Conteúdo |
|---|---|
| [`docs/01_analise_exploratoria.md`](docs/01_analise_exploratoria.md) | fenótipos, UMAP, heatmap, marcadores, interpretação |
| [`docs/08_fenotipos_avancados.md`](docs/08_fenotipos_avancados.md) | missingness, associações, AFP, Kaplan-Meier, Cox, `cox.zph` |
| [`docs/09_auditorias_benchmarks.md`](docs/09_auditorias_benchmarks.md) | auditorias, benchmarks e compliance |
| [`docs/10_grupos_alternativos.md`](docs/10_grupos_alternativos.md) | recidivas, subtipos histológicos, grupos virais e TSS |
| [`docs/11_grupos_separados.md`](docs/11_grupos_separados.md) | análises separadas: recidivas, histologia, vírus e TSS |
| [`docs/02_ComBat.md`](docs/02_ComBat.md) | o que é, modelo, funcionamento, aplicação/limitação na deconvolução |
| [`docs/03_logFC_limma.md`](docs/03_logFC_limma.md) | modelo linear, contrastes, `topTable`, exemplos reais |
| [`docs/04_volcano_plots.md`](docs/04_volcano_plots.md) | construção e interpretação dos 3 volcanos principais + complementares |
| [`docs/05_relatorio_final.md`](docs/05_relatorio_final.md) | síntese, métodos, resultados, limitações, conclusões |
| [`docs/06_dicionario_dados.md`](docs/06_dicionario_dados.md) | dicionário de colunas, fontes e convenções |
| [`docs/07_reproducibilidade.md`](docs/07_reproducibilidade.md) | passo a passo para reproduzir |

### Auditorias e benchmarks

- `auditorias/AUDITORIA_DADOS.md`
- `auditorias/AUDITORIA_VOLCANO_12.md`
- `auditorias/BENCHMARKS.md`
- `auditorias/COMPLIANCE.md`

### Volcano plots principais

| Figura | Comparação |
|---|---|
| `resultados/figuras/Volcano_01_LIHC_TCGA_tumor_vs_adjacente.png` | LIHC (TCGA) |
| `resultados/figuras/Volcano_02_GTEx_masculino_vs_feminino.png` | GTEx |
| `resultados/figuras/Volcano_03_TCGA_tumor_vs_GTEx_normal.png` | TCGA vs. GTEx |

Também há 9 volcanos complementares, todos em PNG.

---

## 6. Como reproduzir

### 6.1 Dependências

```r
# R >= 4.5
install.packages(c("data.table","ggplot2","ggrepel","pheatmap","dplyr",
                   "tidyr","scales","RColorBrewer","matrixStats","uwot",
                   "survival","survminer","broom","patchwork","knitr","png"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("limma","sva","org.Hs.eg.db","AnnotationDbi","edgeR"))
```

Versões testadas em `docs/sessionInfo.txt`.

### 6.2 Dados

Os arquivos pequenos já estão em `dados/raw/`. Para baixar novamente:

```bash
Rscript scripts/00_download_dados.R
```

Para baixar o transcriptoma completo (~1,3 GB):

```bash
BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE Rscript scripts/00_download_dados.R
```

Alternativa: use a matriz processada já incluída
(`dados/processed/expressao_log2tpm_filtrada.tsv.gz`) e pule a etapa 00 de
leitura do TPM bruto. Nesse caso, os scripts 01–05 funcionam diretamente.

### 6.3 Execução

```bash
cd seminários_ryan

# Opção A: pipeline completo (usa TPM em $TEMP/hepato_full/... se existir)
Rscript scripts/run_all.R

# Opção B: passo a passo
TCGA_GTEX_TPM="C:/caminho/TcgaTargetGtex_rsem_gene_tpm.gz" \
  Rscript scripts/00_preparar_dados.R
Rscript scripts/01_analise_exploratoria.R
Rscript scripts/02_expressao_diferencial.R
Rscript scripts/03_volcano_plots.R
Rscript scripts/04_ComBat_demonstracao.R
Rscript scripts/05_fenotipos_avancados.R
Rscript scripts/06_grupos_alternativos.R
Rscript scripts/07_grupos_separados.R
Rscript scripts/10_auditoria_dados_e_volcano.R
Rscript scripts/11_benchmarks.R
Rscript scripts/12_compliance.R
```

Os logs de execução ficam em `resultados/qc/`.

Para gerar novamente o relatório consolidado em Markdown:

```bash
python scripts/09_gerar_relatorio_md.py
```

O arquivo `relatorio_completo.md` é apenas uma conveniência de leitura; a
fonte editável são os arquivos em `docs/`. O projeto não gera HTML.

---

## 7. Limitações

1. **Confundimento TCGA × GTEx:** estudo e condição são perfeitamente
   confundidos; o `logFC` tumor vs. normal combina biologia e lote.
2. **Composição celular:** bulk mede a média de tipos celulares; a
   deconvolução não foi realizada.
3. **TPM:** escala relativa; para contagens, `voom`/modelos binomial negativo
   seriam mais adequados.
4. **Subgrupos pequenos:** estágio IV (n = 5) e grau G4 (n = 11) têm poder
   limitado.
5. **Idade no GTEx:** faixa convertida para ponto médio.
6. **Sem validação experimental:** não houve qPCR, IHC ou coorte independente.
7. **ComBat não resolve confundimento:** a demonstração mostra sobrecorreção.
8. **Missingness clínico:** invasão vascular (70,7%), sorologia viral (55,8%),
   fibrose (43,1%), Child-Pugh (35,8%) e AFP (25,2%) limitam os modelos que
   usam essas variáveis.

---

## 8. Declaração de uso de inteligência artificial

Este material didático foi elaborado com auxílio de ferramentas de IA
generativa para organização do pipeline, revisão de código em R e redação
técnica. Todos os resultados foram executados, conferidos e validados a partir
dos dados públicos descritos; a responsabilidade final pelo conteúdo é do
autor do seminário.

---

## 9. Referências

- Ritchie, M. E. et al. **limma powers differential expression analyses for
  RNA-sequencing and microarray studies.** *Nucleic Acids Research*, 2015.
- Johnson, W. E.; Li, C.; Rabinovic, A. **Adjusting batch effects in microarray
  expression data using empirical Bayes methods.** *Biostatistics*, 2007.
- Leek, J. T. et al. **The sva package for removing batch effects and other
  unwanted variation in high-throughput experiments.** *Bioinformatics*, 2012.
- Zhang, Y.; Parmigiani, G.; Johnson, W. E. **ComBat-seq: batch effect
  adjustment for RNA-seq count data.** *NAR Genomics and Bioinformatics*, 2020.
- Goldman, M. J. et al. **Visualizing and interpreting cancer genomics data via
  the Xena platform.** *Nature Biotechnology*, 2020.
- GTEx Consortium. **The GTEx Consortium atlas of genetic regulatory effects
  across human tissues.** *Science*, 2020.
- Cancer Genome Atlas Research Network. **Comprehensive and integrative
  genomic characterization of hepatocellular carcinoma.** *Cell*, 2017.
- Newman, A. M. et al. **Determining cell type abundance and expression from
  bulk tissues with digital cytometry.** *Nature Biotechnology*, 2019.
- Wang, X. et al. **Bulk tissue cell type deconvolution with multi-subject
  single-cell expression reference.** *Nature Communications*, 2019.
- Korsunsky, I. et al. **Fast, sensitive and accurate integration of
  single-cell data with Harmony.** *Nature Methods*, 2019.
