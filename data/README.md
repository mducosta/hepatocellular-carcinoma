# Dados de entrada

O arquivo de dados **não está versionado** neste repositório (por ser grande e de
terceiros). Para reproduzir a análise, baixe o dataset público do **UCSC Xena**
e salve-o nesta pasta com o nome:

```
liver_lip_aterosclerose.tsv
```

> O script `pipeline_hepato.R` aceita as variantes `.tsv`, `.csv` ou `.xlsx`
> e os nomes `liver_lip_aterosclerose.*` ou `liver_lip_aterosclero.*`.

## Como baixar (UCSC Xena)

1. Acesse <https://xenabrowser.net/>.
2. Adicione os dois conjuntos de dados (datasets):
   - **TCGA Liver Hepatocellular Carcinoma (LIHC)** — gene expression RNAseq;
   - **GTEx Liver** — gene expression RNAseq.
3. Na visualização, selecione **Gene Expression**.
4. Em **File → Download**, baixe o arquivo TSV completo.
5. Renomeie/salve como `data/liver_lip_aterosclerose.tsv`.

## Estrutura esperada do arquivo

O arquivo deve conter colunas de metadados (ex.: `sample`, `sample_type`,
`study`, `TCGA_GTEX_main_category`, `primary disease or tissue`) **e** colunas
de expressão gênica (genes, em escala log2). O pipeline detecta e separa essas
colunas automaticamente.
