# Seminário Ryan — TCGA-LIHC × GTEx (fígado)

Análise exploratória, clínica e de **expressão diferencial** do
**hepatocarcinoma (TCGA-LIHC)** em comparação com **fígado normal (GTEx)**,
com foco em **covariáveis clínicas**, **correção de lote (ComBat)**,
**deconvolução celular** e **robustez do sinal**.

> Pasta: `seminários_ryan/`, dentro do repositório
> `hepatocellular-carcinoma`.
>
> **Artigo final:** [`DOCX`](artigo_final/artigo_expressao_diferencial_covariaveis_deconvolucao.docx)
> (sem resumo, introdução, cabeçalho, rodapé ou paginação).

---

## 1. Pergunta e desenho

**Pergunta:** como o transcriptoma do hepatocarcinoma difere do fígado normal
e do tecido adjacente, e até onde esse sinal sobrevive ao ajuste por sexo,
idade, estágio, grau, sorologia, Child-Pugh, histologia, centro TCGA (TSS),
lote e composição celular?

O desenho combina:

| Bloco | Comparação |
|---|---|
| TCGA × TCGA | Tumor primário vs. tecido adjacente (não pareado e pareado) |
| TCGA × TCGA | Estágio III–IV vs. I–II; grau G3–G4 vs. G1–G2 |
| TCGA × TCGA | Sexo, idade, AFP, Child-Pugh, histologia, fibrose e invasão vascular |
| GTEx × GTEx | Sexo e idade |
| TCGA × GTEx | Tumor primário vs. fígado normal; adjacente vs. fígado normal |
| Via KEGG hsa05417 | Reanálise com covariáveis e composição |
| Sensibilidade | 17 cenários de ComBat (estudo e TSS) |
| Deconvolução | EPIC, quanTIseq e xCell |

---

## 2. Dados

| Grupo | n |
|---|---:|
| TCGA-LIHC — tumor primário | 369 |
| TCGA-LIHC — tecido adjacente | 50 |
| GTEx — fígado normal (após deduplicação) | 109 |
| **Total da análise principal** | **528** |

- **Espécie:** *Homo sapiens*; tecido hepático.
- **Expressão:** UCSC Xena / Toil, `TcgaTargetGtex_rsem_gene_tpm`
  (Toil recompute **v2016-09-03**), **log2(TPM + 0,001)**.
- **Genes após filtro:** 12.401.
- **Painel KEGG hsa05417:** 213 genes (consulta em 17/10/2026).
- **Clínico:** TCGA GDC/`LIHC_clinicalMatrix` (**v2019-12-06**); GTEx
  **Analysis v8** (fenótipos de sexo e faixa etária).
- **Excluídos:** 2 tumores recidivados (`-02`) e 1 segundo fragmento do
  doador GTEX-WYVS (pseudorreplicação). O metadado bruto tem 529 amostras
  (110 GTEx); a análise principal usa 528.
- **Acesso aos dados:** 17 de outubro de 2026.
- **Acesso aos artigos citados:** 26 de outubro de 2026.

> O transcriptoma completo bruto (~1,3 GB) não está versionado. A matriz
> processada (17 MB) está em
> `dados/processed/expressao_log2tpm_filtrada.tsv.gz`.

**TPM** significa *transcripts per million*: medida relativa em que a soma
dos transcritos de uma amostra é normalizada para 1 milhão. Não é contagem
bruta. O valor 0,001 evita logaritmo de zero.

---

## 3. Estrutura da pasta

```text
seminários_ryan/
├── LEIA_ME.md
├── relatorio_consolidado_seminario.md
├── dependencias_R.txt
├── apresentacao_seminario_hepatocarcinoma.pptx
├── scripts/
│   ├── etapa_00_baixar_dados_brutos.R
│   ├── etapa_01_preparar_matrizes_e_metadados.R
│   ├── etapa_02_analise_exploratoria_fenotipos.R
│   ├── etapa_03_expressao_diferencial_limma.R
│   ├── etapa_04_volcano_plots.R
│   ├── etapa_05_correcao_lote_ComBat.R
│   ├── etapa_06_fenotipos_clinicos_avancados.R
│   ├── etapa_07_grupos_alternativos.R
│   ├── etapa_08_analises_grupos_separados.R
│   ├── etapa_09_gerar_relatorio_consolidado.py
│   ├── etapa_10_auditoria_dados_volcano.R
│   ├── etapa_11_benchmarks.R
│   ├── etapa_12_compliance.R
│   ├── etapa_13_expressao_diferencial_covariaveis.R
│   ├── etapa_14_ComBat_testes_multiplos.R
│   ├── etapa_15_deconvolucao_celular.R
│   ├── etapa_16_integracao_deconvolucao_expressao.R
│   ├── etapa_16b_via_LA_covariaveis.R
│   ├── etapa_17_auditoria_analises_avancadas.py
│   ├── etapa_18_benchmarks_analises_avancadas.R
│   ├── etapa_19_compliance_analises_avancadas.py
│   ├── utilitarios_analises_avancadas.R
│   └── executar_todas_as_etapas.R
├── dados/
│   ├── raw/
│   └── processed/
├── resultados/
│   ├── tabelas/
│   │   ├── tabela_01_... a tabela_35_...
│   │   ├── expressao_diferencial/
│   │   └── analises_avancadas/
│   │       ├── tabelas/
│   │       │   ├── expressao_diferencial_covariaveis/
│   │       │   ├── ComBat/
│   │       │   ├── deconvolucao/
│   │       │   ├── integracao_composicao/
│   │       │   └── via_LA/
│   │       ├── figuras/
│   │       ├── qc/
│   │       └── dados_derivados/
│   ├── figuras/
│   └── controle_qualidade/
├── auditorias/
│   ├── auditoria_analises_avancadas.md
│   ├── auditoria_confundimento_associacoes.csv
│   ├── auditoria_fluxo_casos.csv
│   ├── auditoria_atenuacao_sinais.csv
│   ├── benchmarks_analises_avancadas.md
│   ├── compliance_analises_avancadas.md
│   ├── auditoria_dados.md
│   ├── auditoria_volcano_ComBat.md
│   ├── benchmarks_pipeline.md
│   └── compliance_pipeline.md
├── docs/
│   ├── analise_exploratoria_fenotipos.md
│   ├── correcao_lote_ComBat.md
│   ├── calculo_logFC_limma.md
│   ├── volcano_plots.md
│   ├── relatorio_final.md
│   ├── dicionario_dados.md
│   ├── reprodutibilidade.md
│   ├── fenotipos_clinicos_avancados.md
│   ├── auditorias_benchmarks.md
│   ├── grupos_alternativos.md
│   ├── analises_grupos_separados.md
│   ├── sessionInfo_analises_exploratorias.txt
│   ├── sessionInfo_analises_avancadas.txt
│   └── versoes_python_analises_avancadas.txt
├── artigo_final/
│   ├── artigo_expressao_diferencial_covariaveis_deconvolucao.docx
│   └── figuras/
└── artigos/                 # referências completas; não versionada
```

> A pasta `artigos/` contém PDFs/XMLs/HTMLs das referências e é mantida
> apenas localmente. Ela não faz parte do repositório.

---

## 4. Resultados principais

### 4.1 Coorte e qualidade

- **Tumor:** 369 amostras; 120 mulheres/249 homens; idade 59,4 ± 13,5 anos.
- **Adjacente:** 50 amostras; 22 mulheres/28 homens; idade 61,7 ± 16,1 anos.
- **GTEx:** 109 amostras; 32 mulheres/77 homens; idade 53,8 ± 11,1 anos.
- **Estágio (tumor):** I = 169; II = 86; III = 85; IV = 5; ausente = 24.
- **Grau (tumor):** G1 = 55; G2 = 177; G3 = 121; G4 = 11; ausente = 5.
- **Sorologia:** HBV = 59; HCV = 19; HBV+HCV = 85; desconhecida = 206.
- **Child-Pugh corrigido:** A = 215; B-C = 22; ausente = 132.
  A versão anterior tratava 132 ausentes como B-C; isso foi corrigido.
- **Histologia:** HCC = 359; misto = 7; fibrolamelar = 3.
- **AFP:** <20 = 146; 20–399 = 66; ≥400 = 64; ausente = 93.
- **Fibrose:** <5 = 132; ≥5 = 78; ausente = 159.
- **Invasão vascular:** micro = 92; macro = 16; ausente = 261.
- **TSS:** 38 centros; DD concentra 144 tumores (39,0%).

### 4.2 Missingness

- Invasão vascular: **70,7%**
- Sorologia viral: **55,8%**
- Fibrose (Ishak): **43,1%**
- Child-Pugh: **35,8%**
- AFP: **25,2%**
- Estágio: **6,5%**
- Grau: **1,4%**

### 4.3 Expressão diferencial

- **Tumor vs. adjacente:** 4.130 DEGs no basal (3.322 up/808 down);
  3.710 no pareado; 3.740 no clínico; 3.736 com TSS.
- **Tumor vs. GTEx:** 3.391 DEGs (2.251 up/1.140 down); comparação
  estruturalmente confundida com estudo/lote.
- **Estágio:** 46 DEGs na coorte completa; **0** nos mesmos casos completos,
  **0** com TSS e **0** com composição. O centro CC concentra 26/90 tumores
  III–IV e 100% dos casos sem Child-Pugh.
- **Grau:** 172 DEGs no basal; 199 no clínico; 247 com TSS; 337 no modelo só
  com TSS; 153 com ajuste por estágio; **0** após composição.
- **Sexo:** 13 DEGs no GTEx; 175 no tumor; 103 no adjacente. Após composição:
  15 no GTEx, 90 no tumor e 22 no adjacente.
- **AFP ≥400:** 455 DEGs no modelo clínico.
- **Sorologia, Child-Pugh, fibrose e invasão vascular:** **0 DEGs** com
  FDR < 0,05 e |logFC| > 1.
- **Histologia:** 740–1.020 DEGs, mas apenas 10 casos “outros” (5 após
  ajuste completo); resultados **não interpretáveis** como biologia de
  subtipo.
- **Via LA (hsa05417):** 49 DEGs no tumor vs. adjacente basal; 46 no clínico;
  21 após composição. No grau clínico, 6 DEGs.

### 4.4 ComBat

- **TCGA vs. GTEx com mod = ~grupo + sexo + idade:** erro de confundimento.
- **Sem covariável biológica:** 3.391 → 1 DEG; R² do lote em PCs:
  0,170 → 0,0002.
- **mean.only:** 1.947 DEGs; R² do lote 0,119.
- **TSS no grau:** 172 → 290 DEGs; R² do lote 0,075 → 0,015.
- **TSS no estágio:** 46 → 0 DEGs; R² do lote 0,073 → 0,017.
- **TSS no tumor vs. adjacente, sem proteger grupo:** 4.130 → 3.038 DEGs;
  R² do lote 0,074 → 0,012.

### 4.5 Deconvolução celular

- **xCell Hepatocytes:** mediana 0,151 no tumor, 0,225 no adjacente e 0,206
  no GTEx (p = 6,2 × 10⁻⁴⁷).
- **xCell Endothelial:** 0,041 no tumor vs. 0,022 no adjacente e ~0 no GTEx.
- **xCell Macrophages:** 0,049 no adjacente vs. 0,016 no tumor e 0,018 no
  GTEx.
- **EPIC otherCells:** não é hepatócito; correlação negativa com ALB, APOA1,
  TTR e CYP3A4.
- **Concordância entre métodos:** Spearman mediano de 0,01; EPIC e xCell não
  são intercambiáveis.
- **Cox ajustado:** xCell Hepatocytes HR = 0,77 por SD (0,65–0,92;
  FDR = 0,029); EPIC otherCells HR = 1,39 (1,13–1,70; FDR = 0,029), mas esse
  estimador não passou no QC.

### 4.6 Sexo e hepatócitos

O sexo altera o transcriptoma hepático (XIST e genes do Y funcionam como
controles positivos). A diferença de hepatócitos por sexo foi pequena:
no tumor, mediana 0,145 em mulheres vs. 0,158 em homens (p = 0,010;
FDR = 0,129); no GTEx, 0,203 vs. 0,207 (p = 0,717); no adjacente, 0,229 vs.
0,223 (p = 0,0016, diferença de pequena magnitude). O ajuste por composição
reduziu, mas não zerou, os DEGs de sexo no tumor. Isso é compatível com
dimorfismo sexual intrínseco ao hepatócito e a outros tipos celulares, como
descrito na literatura (Kočar et al., 2026; DOI 10.1186/s13293-025-00811-7).

---

## 5. Auditoria crítica

Os principais achados de auditoria estão em
[`auditorias/auditoria_analises_avancadas.md`](auditorias/auditoria_analises_avancadas.md):

1. **TCGA × GTEx:** grupo e estudo têm Cramér's V = 1,0; não é efeito puro
   de tumor.
2. **Child-Pugh:** 132 ausentes eram tratados como B-C; B-C cai de 154 para
   22 após a correção.
3. **Estágio × TSS:** o centro CC tem 81,2% de estágio III–IV e 100% sem
   Child-Pugh; a assinatura de estágio não é identificável separadamente do
   centro.
4. **Estágio frágil:** 46 DEGs → 0 nos mesmos casos completos, com TSS ou
   com composição.
5. **Grau:** 172–337 DEGs em modelos diferentes, mas 0 após composição.
6. **Tumor vs. adjacente:** 4.130 → 3.710 pareado → 3.740 clínico; o eixo
   mais robusto.
7. **Deconvolução:** EPIC otherCells não é hepatócito; xCell é escore, não
   fração; a concordância entre métodos é baixa.
8. **Sexo:** 13 DEGs no GTEx e 175 no tumor; 90 no tumor após composição.

---

## 6. Como reproduzir

### 6.1 Dependências

```r
# R >= 4.6
install.packages(c("data.table","ggplot2","ggrepel","pheatmap","dplyr",
                   "tidyr","scales","RColorBrewer","matrixStats","uwot",
                   "survival","survminer","broom","patchwork","knitr","png",
                   "vegan","car"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("limma","sva","edgeR","EPIC","quantiseqr","xCell",
                       "GSVA","org.Hs.eg.db","AnnotationDbi"))
```

Versões testadas em `docs/sessionInfo_analises_avancadas.txt` e
`docs/versoes_python_analises_avancadas.txt`.

### 6.2 Dados

Os arquivos pequenos já estão em `dados/raw/`. Para baixar novamente:

```bash
Rscript scripts/etapa_00_baixar_dados_brutos.R
```

Para baixar o transcriptoma completo (~1,3 GB):

```bash
BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE Rscript scripts/etapa_00_baixar_dados_brutos.R
```

Alternativa: use a matriz processada já incluída
(`dados/processed/expressao_log2tpm_filtrada.tsv.gz`).

### 6.3 Execução

```bash
cd seminários_ryan

# Pipeline completo, incluindo análises avançadas, figuras e artigo
Rscript scripts/executar_todas_as_etapas.R

# Ou passo a passo (exemplo)
Rscript scripts/etapa_13_expressao_diferencial_covariaveis.R
Rscript scripts/etapa_14_ComBat_testes_multiplos.R
Rscript scripts/etapa_15_deconvolucao_celular.R
Rscript scripts/etapa_16_integracao_deconvolucao_expressao.R
Rscript scripts/etapa_16b_via_LA_covariaveis.R
python scripts/etapa_17_auditoria_analises_avancadas.py
Rscript scripts/etapa_18_benchmarks_analises_avancadas.R
python scripts/etapa_19_compliance_analises_avancadas.py
```

Os logs ficam em `resultados/controle_qualidade/` e
`resultados/analises_avancadas/qc/`.

### 6.4 Artigos citados

Os PDFs/XMLs/HTMLs das referências não são versionados. A pasta `artigos/` é
mantida apenas na cópia local. O artigo de Kočar et al. (2026), usado para
justificar o dimorfismo sexual em hepatócitos, está entre as referências
citadas no artigo final.

---

## 7. Limitações

1. **Confundimento TCGA × GTEx:** estudo e condição são perfeitamente
   confundidos; o logFC combina biologia e lote.
2. **Composição celular:** bulk mede a média de tipos celulares; a
   deconvolução é inferência, não medida direta.
3. **TPM:** escala relativa; não usamos contagens brutas, voom ou modelos
   binomial negativo.
4. **Missingness:** invasão vascular 70,7%, sorologia 55,8%, fibrose 43,1%,
   Child-Pugh 35,8% e AFP 25,2%; casos completos não são representativos.
5. **Subgrupos pequenos:** estágio IV (n = 5), G4 (n = 11) e histologia
   “outros” (n = 10).
6. **TSS:** proxy imperfeito de lote; confundido com estágio, grau e
   histologia.
7. **Deconvolução:** métodos discordantes; EPIC otherCells não é hepatócito;
   xCell é escore, não fração.
8. **Idade no GTEx:** faixa convertida ao ponto médio.
9. **Sem validação experimental:** não houve qPCR, IHC, Western blot ou
   coorte independente.
10. **Sexo:** sem dados de menopausa, hormônios, genótipo do HCV, carga viral
    ou tratamento; a análise é incompleta.
11. **Dados observacionais:** associações não implicam causalidade.

---

## 8. Documentação e auditorias

| Documento | Conteúdo |
|---|---|
| [`docs/analise_exploratoria_fenotipos.md`](docs/analise_exploratoria_fenotipos.md) | fenótipos, UMAP, heatmap, marcadores |
| [`docs/correcao_lote_ComBat.md`](docs/correcao_lote_ComBat.md) | ComBat, modelo e limitações |
| [`docs/calculo_logFC_limma.md`](docs/calculo_logFC_limma.md) | logFC, contrastes e `topTable` |
| [`docs/volcano_plots.md`](docs/volcano_plots.md) | construção e leitura dos volcanos |
| [`docs/relatorio_final.md`](docs/relatorio_final.md) | síntese do seminário |
| [`docs/dicionario_dados.md`](docs/dicionario_dados.md) | dicionário de colunas |
| [`docs/reprodutibilidade.md`](docs/reprodutibilidade.md) | passo a passo |
| [`docs/fenotipos_clinicos_avancados.md`](docs/fenotipos_clinicos_avancados.md) | missingness, Cox, Kaplan-Meier |
| [`docs/grupos_alternativos.md`](docs/grupos_alternativos.md) | recidivas, histologia, vírus, TSS |
| [`docs/analises_grupos_separados.md`](docs/analises_grupos_separados.md) | análises separadas |
| [`auditorias/auditoria_analises_avancadas.md`](auditorias/auditoria_analises_avancadas.md) | auditoria crítica completa |
| [`auditorias/benchmarks_analises_avancadas.md`](auditorias/benchmarks_analises_avancadas.md) | tempos, memória e concordância |
| [`auditorias/compliance_analises_avancadas.md`](auditorias/compliance_analises_avancadas.md) | checklist de conformidade |

---

## 9. Declaração de uso de inteligência artificial

Em conformidade com a Portaria CNPq n° 2.664/2026, declara-se que ferramentas
de IA generativa foram utilizadas como assistentes de programação, organização
do pipeline, revisão de código e redação técnica. Todas as rotinas foram
executadas, conferidas e validadas contra os arquivos de saída; a interpretação
dos resultados, a definição do desenho analítico e a responsabilidade
intelectual pelo conteúdo permanecem sob responsabilidade do autor. Nenhum
dado identificável de paciente foi usado.

---

## 10. Referências principais

- Kočar E, Blagotinšek Cokan K, Kreft T, Režen T, Rozman D. The importance of
  sex dimorphism in liver metabolism and progressive liver diseases. Biology
  of Sex Differences. 2026;17(1):14. doi:10.1186/s13293-025-00811-7.
- Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression
  analyses for RNA-sequencing and microarray studies. Nucleic Acids Research.
  2015;43(7):e47. doi:10.1093/nar/gkv007.
- Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. The sva package for
  removing batch effects and other unwanted variation in high-throughput
  experiments. Bioinformatics. 2012;28(6):882-883.
  doi:10.1093/bioinformatics/bts034.
- Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular
  heterogeneity landscape. Genome Biology. 2017;18(1):220.
  doi:10.1186/s13059-017-1349-1.
- Wang X, Park J, Susztak K, Zhang NR, Li M. Bulk tissue cell type
  deconvolution with multi-subject single-cell expression reference. Nature
  Communications. 2019;10:380. doi:10.1038/s41467-018-08023-x.
- Finotello F, Mayer C, Plattner C, et al. Molecular and pharmacological
  modulators of the tumor immune contexture revealed by deconvolution of
  RNA-seq data. Genome Medicine. 2019;11(1):34.
  doi:10.1186/s13073-019-0638-6.

A lista completa de referências está acima e no artigo final. Os arquivos
integrais (PDFs/XMLs/HTMLs) não são versionados no repositório.
