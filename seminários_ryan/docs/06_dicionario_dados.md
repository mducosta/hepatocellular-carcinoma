# 6. Dicionário de dados e arquivos

> Arquivo: `docs/06_dicionario_dados.md`.

---

## 6.1 Fontes dos dados

| Arquivo | Fonte | URL |
|---|---|---|
| `dados/raw/LIHC_clinicalMatrix` | UCSC Xena / TCGA LIHC clinical matrix | `https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/LIHC_clinicalMatrix` |
| `dados/raw/TCGA-LIHC.clinical.tsv.gz` | UCSC Xena / GDC TCGA-LIHC | `https://gdc.xenahubs.net/download/TCGA-LIHC.clinical.tsv.gz` |
| `dados/raw/GTEX_phenotype.gz` | UCSC Xena / GTEx phenotype | `https://toil.xenahubs.net/download/GTEX_phenotype.gz` |
| `dados/raw/TcgaTargetGTEX_phenotype.txt.gz` | UCSC Xena / Toil phenotype | `https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz` |
| `dados/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` | GTEx Portal v8 | `https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` |
| `dados/raw/liver.tsv` | painel hepático (213 genes da via KEGG hsa05417 + metadados); cópia local do projeto principal | repositório `hepatocellular-carcinoma` |
| `TcgaTargetGtex_rsem_gene_tpm.gz` | UCSC Xena / Toil (transcriptoma completo) | `https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz` |

O transcriptoma completo (≈1,3 GB) **não** foi versionado dentro de
`seminários_ryan/`; ele foi processado e a matriz filtrada (17 MB) foi salva em
`dados/processed/expressao_log2tpm_filtrada.tsv.gz`. Para reproduzir do zero,
use `scripts/00_download_dados.R` ou baixe manualmente pela URL acima e defina
a variável de ambiente `TCGA_GTEX_TPM` com o caminho do arquivo.

---

## 6.2 Seleção das amostras

- **TCGA-LIHC tumor primário:** `_sample_type == "Primary Tumor"` e
  `_study == "TCGA"` → 369 amostras.
- **TCGA-LIHC tecido adjacente:** `_sample_type == "Solid Tissue Normal"` →
  50 amostras.
- **GTEx fígado normal:** `_study == "GTEX"`, `_primary_site == "Liver"`,
  `_sample_type == "Normal Tissue"` → 110 amostras.
- **Recidivas:** 2 amostras `Recurrent Tumor` (`TCGA-DD-AACA-02` e
  `TCGA-ZS-A9CF-02`) foram **excluídas** da análise principal.
- **Total:** 531 amostras no metadado completo; 529 na análise principal.

---

## 6.3 `metadata_amostras.tsv` e `metadata_analise.tsv`

| Coluna | Descrição |
|---|---|
| `sample` | identificador da amostra (TCGA barcode ou GTEx sample ID) |
| `samples` | cópia do identificador (coluna original do liver.tsv) |
| `detailed_category` | categoria detalhada do Xena (`Liver Hepatocellular Carcinoma`, etc.) |
| `DFI`, `DFI.time` | doença livre de intervalo e tempo (quando disponível) |
| `DSS`, `DSS.time` | sobrevida específica da doença e tempo |
| `OS`, `OS.time` | sobrevida global e tempo (dias) |
| `PFI`, `PFI.time` | intervalo livre de progressão e tempo |
| `sex_liver` | sexo original do liver.tsv (`_gender`) |
| `sample_type` | tipo de amostra (`Primary Tumor`, `Solid Tissue Normal`, `Normal Tissue`, `Recurrent Tumor`) |
| `study` | `TCGA` ou `GTEX` |
| `main_category` | categoria principal do Xena (`TCGA Liver Hepatocellular Carcinoma` ou `GTEX Liver`) |
| `condition` | `Tumor`, `Adjacent`, `Normal` ou `Recurrent` |
| `group` | `TCGA_Tumor`, `TCGA_Adjacent`, `GTEx_Normal`, `TCGA_Recurrent` |
| `dataset` | `TCGA-LIHC` ou `GTEx-Liver` |
| `subject` | identificador do paciente/doador (TCGA: `TCGA-XX-YYYY`; GTEx: `GTEX-XXXX`) |
| `sex_tcga` | sexo do clínico TCGA |
| `age` | idade em anos (TCGA: idade ao diagnóstico; GTEx: ponto médio da faixa) |
| `stage_raw` | estágio AJCC original (texto) |
| `stage` | estágio harmonizado (`I`, `II`, `III`, `IV`) |
| `grade_raw` | grau histológico original (`G1`–`G4`, `Not Reported`) |
| `grade` | grau harmonizado (`G1`–`G4`) |
| `histology` | tipo histológico |
| `vital_status` | status vital |
| `days_to_death` | dias até o óbito |
| `days_to_last_followup` | dias de seguimento |
| `fibrosis_ishak` | escore de fibrose de Ishak |
| `afp` | alfa-fetoproteína (valor laboratorial) |
| `viral_hepatitis` | sorologia para hepatites virais |
| `child_pugh` | classificação de Child-Pugh |
| `sex_gtex` | sexo do GTEx (1 = masculino, 2 = feminino) |
| `age_gtex` | idade do GTEx (ponto médio da faixa) |
| `sex` | sexo harmonizado (`female`/`male`) |
| `age_group` | `<65` ou `>=65` |
| `age_group60` | `<60` ou `>=60` |
| `stage_group` | `I-II` ou `III-IV` |
| `grade_group` | `G1-G2` ou `G3-G4` |
| `tss` | código do centro TCGA (Tissue Source Site) |
| `viral_status` | `None/Other`, `HBV`, `HCV`, `HBV+HCV` ou `Unknown` |
| `afp_log2` | log2(AFP + 1) |
| `afp_high` | AFP `<400` ou `>=400` ng/mL |
| `afp_group` | AFP `<20`, `20-399` ou `>=400` |
| `fibrosis_group` | Ishak `0`, `1-4`, `5` ou `6` |
| `fibrosis_advanced` | Ishak `<5` ou `>=5` |
| `child_pugh_group` | Child-Pugh `A` ou `B-C` |
| `histology_group` | `HCC` ou `Other` |
| `vascular_invasion` | `None/Unknown`, `Micro` ou `Macro` |
| `race`, `race_group` | raça declarada e agrupada |
| `ethnicity`, `ethnicity_group` | etnia declarada e agrupada |
| `pathologic_t`, `pathologic_n`, `pathologic_m` | T/N/M patológicos |
| `ajcc_t`, `ajcc_n`, `ajcc_m` | T/N/M AJCC (GDC) |
| `residual_tumor` | R0, R1, R2 ou RX |
| `neoplasm_status` | status neoplásico no seguimento |

> Observação: para as 50 amostras adjacentes, `stage` e `grade` refletem o
> tumor do mesmo paciente, não o tecido adjacente em si. As análises de
> estágio/grau usam apenas `group == "TCGA_Tumor"`.

---

## 6.4 Matrizes de expressão

### `expressao_log2tpm_filtrada.tsv.gz`

- Formato: TSV gzipado, primeira coluna `gene` (símbolo HUGO), demais colunas
  amostras.
- Valores: log2(TPM + 0,001).
- Dimensões: 12.401 genes × 529 amostras.
- Filtro: expressão ≥ 1 em ≥ 20% das amostras; variância > 0.

### `expressao_LA_pathway.tsv.gz`

- Formato: TSV gzipado, primeira coluna `gene`.
- Genes: 213 genes da via KEGG hsa05417 (Lipídios e Aterosclerose).
- Valores: log2(TPM + 0,001).
- Amostras: 531 (inclui recidivas, mas as análises principais usam 529).

### `resumo_preparacao.txt`

- Resumo textual da etapa 00, com contagens, dimensões e sessão R.

---

## 6.5 Tabelas de expressão diferencial

Arquivos em `resultados/tabelas/DE/`:

| Arquivo | Conteúdo |
|---|---|
| `DEG_<comparacao>_completo.csv` | todos os genes testados, com estatísticas completas |
| `DEG_<comparacao>_significativos.csv` | apenas DEGs (FDR < 0,05 e \|logFC\| > 1) |
| `TOP20_<comparacao>.csv` | top 20 up + top 20 down |
| `RESUMO_DEGs_TODAS_COMPARACOES.csv` | contagens de up/down por comparação |
| `COMPARACAO_ComBat_antes_depois.csv` | comparação dos DEGs antes/depois do ComBat |

### Colunas das tabelas de DEG

| Coluna | Descrição |
|---|---|
| `gene` | símbolo HUGO |
| `comparison` | nome da comparação |
| `direction` | direção do contraste (ex.: `TCGA_Tumor > GTEx_Normal`) |
| `logFC` | log2 fold change |
| `AveExpr` | expressão média do gene (log2) |
| `t` | estatística t moderada do limma |
| `P.Value` | p-valor nominal |
| `adj.P.Val` | p-valor ajustado por BH (FDR) |
| `B` | log-odds bayesiano de diferença |
| `regulation` | `Up`, `Down` ou `NS` |

### Nomes das comparações

- `TCGA_Tumor_vs_TCGA_Adjacent`
- `TCGA_Tumor_vs_Adjacent_PAREADO`
- `TCGA_Tumor_vs_GTEx_Normal`
- `TCGA_Adjacent_vs_GTEx_Normal`
- `GTEx_Male_vs_Female`
- `GTEx_Age60plus_vs_lt60`
- `TCGA_Stage_III_IV_vs_I_II`
- `TCGA_Grade_G3_G4_vs_G1_G2`
- `TCGA_Male_vs_Female`
- `TCGA_Age65plus_vs_lt65`
- `LApathway_TCGA_Tumor_vs_GTEx_Normal`
- `LApathway_TCGA_Tumor_vs_TCGA_Adjacent`
- `TCGA_Tumor_vs_GTEx_Normal_ComBat` (demonstração de sobrecorreção)

---

## 6.6 Tabelas clínicas avançadas

| Arquivo | Conteúdo |
|---|---|
| `13_fenotipos_missingness.csv` | missingness das variáveis clínicas |
| `14_fenotipos_testes.csv` | testes de associação e Cramér's V |
| `15_AFP_por_subgrupo.csv` | resumo de AFP por estágio, grau, sexo e vírus |
| `16_KM_logrank.csv` | log-rank das curvas de sobrevida |
| `17_cox_univariado.csv` | Cox univariado |
| `18_cox_multivariado.csv` | Cox multivariado (Modelos A e B) |
| `19_correlacao_clinica.csv` | correlações de Spearman |
| `20_cox_ph_assumption.csv` | `cox.zph` (riscos proporcionais) |

## 6.7 Figuras

| Arquivo | Conteúdo |
|---|---|
| `01_amostras_por_grupo.png` | distribuição das amostras |
| `02_sexo_por_grupo.png` | composição por sexo |
| `03_idade_por_grupo.png` | idade por grupo |
| `04_estagio_TCGA.png` | estágio AJCC |
| `05_grau_TCGA.png` | grau histológico |
| `06_estagio_x_grau_TCGA.png` | estágio × grau |
| `07_UMAP_grupo.png` | UMAP por grupo |
| `08_UMAP_estudo.png` | UMAP por estudo |
| `09_UMAP_sexo_TCGA.png` | UMAP por sexo (TCGA tumor) |
| `10_UMAP_estagio_grau_TCGA.png` | UMAP por estágio/grau |
| `11_heatmap_top50_variaveis.png` | heatmap top 50 genes variáveis |
| `12_boxplot_marcadores.png` | boxplots de marcadores |
| `13_UMAP_ComBat_antes_depois.png` | UMAP antes/depois do ComBat |
| `14_AFP_ComBat_antes_depois.png` | AFP antes/depois do ComBat |
| `15_AFP_por_estagio_grau.png` | AFP por estágio e grau |
| `16_missingness_clinico.png` | missingness das variáveis clínicas |
| `17_KM_OS_estagio.png` | Kaplan-Meier por estágio |
| `18_KM_OS_grau.png` | Kaplan-Meier por grau |
| `19_KM_OS_sexo.png` | Kaplan-Meier por sexo |
| `20_KM_OS_idade.png` | Kaplan-Meier por idade |
| `21_KM_OS_AFP.png` | Kaplan-Meier por AFP |
| `22_KM_OS_fibrose.png` | Kaplan-Meier por fibrose |
| `23_KM_OS_painel.png` | painel de curvas KM |
| `24_forest_cox_multivariado.png` | forest plot do Cox multivariado |
| `25_correlacao_clinica.png` | heatmap de correlações clínicas |
| `Volcano_01_LIHC_TCGA_tumor_vs_adjacente.png` | volcano principal LIHC |
| `Volcano_02_GTEx_masculino_vs_feminino.png` | volcano principal GTEx |
| `Volcano_03_TCGA_tumor_vs_GTEx_normal.png` | volcano principal TCGA vs. GTEx |
| `Volcano_04...12...` | volcanos complementares |

---

## 6.8 Convenções

- **logFC:** sempre em log2; positivo = maior no grupo do numerador do
  contraste.
- **FDR:** Benjamini-Hochberg.
- **Significância:** FDR < 0,05 e |logFC| > 1, salvo indicação contrária.
- **Idade:** anos; no GTEx, ponto médio da faixa.
- **Estágio:** AJCC harmonizado em I–IV.
- **Grau:** G1–G4.
- **Cores dos grupos:** TCGA tumor (vermelho), TCGA adjacente (laranja),
  GTEx normal (azul).
