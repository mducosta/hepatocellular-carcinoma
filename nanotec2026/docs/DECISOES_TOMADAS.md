# Decisões tomadas — Nanotec 2026 (GINAB/UENF)

Estudo: **Evidência transcriptômica de remodelamento da matriz extracelular (MEC)
e superexpressão de MMP9/MMPs/TIMPs no hepatocarcinoma (HCC)**.

Toda decisão abaixo foi tomada com base em **evidência observada nos próprios
dados** ou em justificativa metodológica explícita. Nenhum dado foi inventado.

---

## 1. Fonte e download

| Decisão | Justificativa (evidência) |
|---|---|
| Hub **Toil** (`https://toil.xenahubs.net`) | Recomendado no protocolo; é a fonte harmonizada TCGA+GTEx usada em estudos de expressão. |
| Dataset **`TcgaTargetGtex_rsem_gene_tpm`** | Recomendado no protocolo (TPM harmonizado). Confirmamos via `UCSCXenaTools::XenaQuery` que a URL canônica é `.../TcgaTargetGtex_rsem_gene_tpm.gz`. |
| Fenótipo **`TcgaTargetGTEX_phenotype.txt`** | Necessário para filtrar tecido/tipo de amostra. |
| Sobrevivência **`TCGA_survival_data`** | Contém OS/OS.time (e DSS/PFI), usado no KM/Cox. |
| **Baixar o transcriptoma completo (1,32 GB) em vez de só o painel** | Dois motivos: (i) os endpoints de *single-gene* do Toil (`.../download/<dataset>/<gene>`) retornam **HTTP 403 (AccessDenied)** no S3 — só os arquivos físicos `.gz` estão acessíveis; (ii) **GSVA/ssGSEA exige o ranking de TODOS os genes por amostra** como fundo. |

### Evidência de conectividade (registrada durante a auditoria)
- `.../TcgaTargetGtex_rsem_gene_tpm.gz` → **HTTP 206** (acessível, 1.323.254.426 bytes).
- `.../TcgaTargetGTEX_phenotype.txt.gz` → **HTTP 206** (135.753 bytes).
- `.../TCGA_survival_data` → **HTTP 206** (406.381 bytes).
- `.../download/TcgaTargetGtex_rsem_gene_tpm/MMP9` (single gene) → **HTTP 403**.

---

## 2. Descoberta crítica: IDs Ensembl e mapeamento

O arquivo `TcgaTargetGtex_rsem_gene_tpm` **não usa símbolos HUGO**; as linhas são
Ensembl gene IDs com versão (ex.: `ENSG00000100985.16`).

**Decisão:** mapear Ensembl → símbolo localmente com `org.Hs.eg.db` (sem rede
externa), removendo o sufixo de versão e agregando isoformas/versões pela média.

**Evidência:** 60.498 linhas Ensembl lidas; 35.545 IDs mapeados para 36.301
símbolos únicos; 24.953 linhas sem símbolo (lncRNAs/pseudogenes) descartadas.
Todos os **41 genes do painel** foram recuperados após o mapeamento.

---

## 3. Filtro de amostras de fígado

**Evidência do fenótipo** (cruzamento `_study` × `_sample_type` para `_primary_site == "Liver"`):

| Grupo | Critério | n |
|---|---|---|
| **Tumor** (LIHC) | `_study=TCGA` & `_sample_type=Primary Tumor` | 369 |
| **Normal** (GTEx) | `_study=GTEX` & `_sample_type=Normal Tissue` | 110 |
| **Adjacent** | `_study=TCGA` & `_sample_type=Solid Tissue Normal` | 50 |
| *Recurrent Tumor* | `_study=TCGA` & `_sample_type=Recurrent Tumor` | 2 → **excluídos** |

**Decisão:** excluir os 2 "Recurrent Tumor" — são recidivas (categoria distinta
de tumor primário); mantê-los misturaria desenho. Total analisado: **529 amostras**.

As 529 amostras foram todas localizadas nas colunas da matriz TPM (interseção 100%).

---

## 4. Escala dos dados

A auditoria encontrou `min = -9.9658`, `max = 14.15`, mediana 3.21, 0% valores > 100.

**Conclusão:** os dados **já estão em log2(TPM + 0.001)** (o piso −9,966 =
log2(0.001)). Portanto **NÃO reaplicamos log2**; usamos limma direto sobre a
escala log2 (adequado para dados contínuos). Isso difere do esperado "TPM linear"
e foi detectado empiricamente.

---

## 5. Modelagem estatística

| Decisão | Justificativa |
|---|---|
| **limma** (lmFit + contrasts.fit + eBayes robusto) | Padrão-ouro para DE em matrizes de expressão contínua (log2); robusto a outliers (`robust=TRUE`). |
| Contraste principal **Tumor vs Normal** | Objetivo central do estudo. |
| Contraste **Tumor vs Adjacent** | Validação *sem* o confundimento de plataforma TCGA vs GTEx (ambos são TCGA). |
| Critério DEG: **FDR < 0,05 e \|log2FC\| > 1** | Convenção; documentada. |
| Correlação MMP9: **Spearman** (não paramétrica) | Evita pressuposto de linearidade/normalidade. |
| Kaplan-Meier: **split pela mediana** da expressão no grupo tumoral | Parcimonioso e transparente. |
| Cox **univariado** (expressão contínua) | O arquivo de sobrevivência não traz covariáveis clínicas (idade/estádio); ajuste multivariado impossível → documentado como limitação. |
| **GSVA com ssGSEA** | "GSVA ou similar"; ssGSEA é rank-based e adequado a um transcriptoma de 36 mil genes. |

---

## 6. Batch effect TCGA vs GTEx

A PCA do painel foi gerada por grupo e por estudo. O contraste **Tumor (TCGA) vs
Normal (GTEx)** tem confundimento de origem (plataforma). Por isso o contraste
**Tumor vs Adjacent** (ambos TCGA) é reportado como verificação: ele confirmou a
mesma direção dos MMPs-chave, indicando que o sinal não é artefato de plataforma.

---

## 7. O que foi descartado e por quê

| Análise | Status | Motivo |
|---|---|---|
| Single-gene download | Descartado | Endpoint retorna HTTP 403 (S3). |
| Cox multivariado | Descartado | Ausência de covariáveis clínicas no `TCGA_survival_data`. |
| Validação proteica | Não realizada | Fora do escopo dos dados baixados (transcriptoma bulk). |
| Nanocarreadores experimentais | **Não incluído** | Sem dados experimentais próprios; apenas perspectivas translacionais. |

---

## 8. Ambiente de execução

- R **4.6.1** (Windows), `Rscript` em `C:/Program Files/R/R-4.6.1/bin`.
- Pacotes: UCSCXenaTools 1.7.0, limma, GSVA 2.6.0, msigdbr 26.1.0, org.Hs.eg.db,
  survival, survminer, ggplot2, ggrepel, pheatmap, data.table.
- **Nota técnica:** o R 4.6.1 (muito recente) apresentou `segfault` intermitente
  ao rodar código via `Rscript -e` no bash (escaping), e o shell do Windows
  (`cmd.exe`) interpreta `||` como operador de shell. Mitigações: (i) todos os
  scripts são **arquivos .R** (nunca `-e`); (ii) o filtro de genes usa **múltiplas
  regras awk** sem `||`; (iii) pipes no Windows usam `shell()` em vez de `system()`.
