# Dados brutos

Arquivos pequenos usados para construir os metadados clínicos:

| Arquivo | Origem | Conteúdo |
|---|---|---|
| `LIHC_clinicalMatrix` | UCSC Xena / TCGA | clínico TCGA-LIHC (idade, sexo, estágio, grau, AFP, fibrose, etc.) |
| `TCGA-LIHC.clinical.tsv.gz` | UCSC Xena / GDC | clínico TCGA-LIHC em formato GDC (validação) |
| `GTEX_phenotype.gz` | UCSC Xena / GTEx | tecido, sexo e amostra GTEx |
| `TcgaTargetGTEX_phenotype.txt.gz` | UCSC Xena / Toil | fenótipo combinado TCGA + GTEx usado na seleção das amostras |
| `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` | GTEx Portal v8 | sexo e faixa etária dos doadores GTEx |
| `liver.tsv` | projeto `hepatocellular-carcinoma` | painel hepático (213 genes da via KEGG hsa05417) + metadados |

## Transcriptoma completo

O arquivo `TcgaTargetGtex_rsem_gene_tpm.gz` (~1,3 GB) **não** está incluído
nesta pasta. Ele pode ser baixado com:

```bash
Rscript ../scripts/00_download_dados.R
# ou
curl -L -o TcgaTargetGtex_rsem_gene_tpm.gz \
  https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz
```

A matriz já processada e filtrada está em
`../processed/expressao_log2tpm_filtrada.tsv.gz`.

## Acesso e licenças

- TCGA: dados públicos do National Cancer Institute (NCI/GDC), sujeitos às
  políticas de uso do TCGA.
- GTEx: dados públicos do GTEx Consortium (dbGaP/GTEx Portal), sujeitos às
  políticas de uso do GTEx.
- UCSC Xena/Toil: camada de distribuição pública dos dados processados.

Consulte os termos originais de cada projeto antes de redistribuir os dados.
