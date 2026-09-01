# Tutorial — Baixar dados no UCSC Xena

Guia passo a passo para baixar o **transcriptoma completo** e os **dados
adicionais** (metilação, CNV, mutações, microRNA, proteômica, scRNA-seq) que
podem ampliar a análise.

> Tudo pode ser feito **somente computacionalmente** (não há etapa experimental
> obrigatória). Os dados são públicos no UCSC Xena / GDC / CPTAC / GEO.

---

## 1. Como baixar o transcriptoma COMPLETO (~20 mil genes)

O arquivo atual (`liver.tsv`) contém apenas os **212 genes da via hsa05417**.
Para análises mais poderosas (GSEA/GSVA/ORA com poder estatístico), baixe o
transcriptoma completo:

1. Acesse <https://xenabrowser.net/>.
2. No menu **"Dataset"**, adicione os conjuntos de **gene expression RNAseq**:
   - **TCGA Liver Cancer (LIHC)** → *gene expression RNAseq (HTSeq - FPKM)*
     (ou *TOIL RSEM*);
   - **GTEx Liver** → *gene expression RNAseq*.
3. A visualização mostrará **~20.000 genes** (transcriptoma inteiro), não só a via.
4. **File → Download** e salve como `liver_full.tsv`.
5. Ajuste o script para apontar para o novo arquivo (ou renomeie para `liver.tsv`).

> ⚠️ O arquivo completo é grande (~100–300 MB). O `.gitignore` ignora os dados
> brutos, **exceto** `dados/raw/liver.tsv` (que é versionado no repositório).

### Bookmark oficial (via hsa05417, usado no estudo)
```
https://xenabrowser.net/?bookmark=337fe0532808c6fc66cf017f13885c4a
```
Para o transcriptoma completo, crie um novo bookmark após adicionar os datasets
completos (File → Save bookmark).

---

## 2. Dados adicionais (todos no Xena/GDC/CPTAC)

| Dado extra | Onde baixar | Análise |
|-----------|-------------|---------|
| **Metilação (Illumina 450k)** | Xena → *TCGA LIHC DNA methylation (450k)* | Epigenética dos DEGs |
| **CNV (GISTIC)** | Xena → *TCGA LIHC copy number (GISTIC2)* | Alterações de cópia dos DEGs |
| **Mutações (MuTect/VarScan)** | Xena → *TCGA LIHC somatic mutation* | Mutações nos genes da via |
| **microRNA (miR)** | Xena → *TCGA LIHC miRNA mature strand* | Regulação pós-transcricional |
| **Proteômica (CPTAC)** | <https://proteomics.cancer.gov> (CPTAC-LIHC) | Validação proteica |
| **scRNA-seq de CHC** | GEO (ex.: GSE149614, GSE125449) | Expressão célula a célula |
| **Coortes de imunoterapia** | IMbrave150 (atezo+bev) — repositórios públicos | Resposta ao tratamento |

### Como baixar metilação/CNV/mutação/miRNA no Xena
1. Em <https://xenabrowser.net/>, digite "LIHC" na busca.
2. Selecione o **tipo de dado** desejado (methylation, copy number, mutation, miRNA).
3. Escolha **"Liver Hepatocellular Carcinoma"** e o tipo de amostra (Primary Tumor).
4. **File → Download** (TSV).
5. Salve na pasta `dados/` e documente o bookmark correspondente.

---

## 3. Coortes independentes para validação

| Coorte | Acesso | Uso |
|--------|--------|-----|
| **ICGC** (LIRI-JP) | <https://dcc.icgc.org> | Validação externa |
| **GEO** | <https://www.ncbi.nlm.nih.gov/geo/> | Buscar "hepatocellular carcinoma expression" |
| **CPTAC-LIHC** | <https://proteomics.cancer.gov> | Validação proteica |

---

## 4. Pipeline após baixar os dados

1. Salve o arquivo na raiz (ex.: `liver_full.tsv`).
2. Rode o pipeline:
   ```bash
   Rscript pipeline_hepato.R
   ```
3. Para os dados adicionais (metilação/CNV/miRNA), crie scripts específicos ou
   use os pacotes `minfi` (metilação), `GISTIC`/`maftools` (CNV/mutação),
   `limma` (miRNA).
