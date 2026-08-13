# Dados de entrada

O arquivo de dados **não está versionado** neste repositório (dados públicos de
terceiros). Para reproduzir a análise, baixe o dataset do **UCSC Xena** e
salve-o na **raiz do repositório** com o nome:

```
liver.tsv
```

> O script `pipeline_hepato.R` aceita também as variantes `.tsv`, `.csv` ou
> `.xlsx` e os nomes `liver.*` ou `liver_lip_aterosclerose.*`.

## Bookmark oficial do dataset (UCSC Xena)

Use o bookmark abaixo para recarregar exatamente o dataset usado neste estudo
(TCGA-LIHC + GTEx Liver, expressão gênica):

```
https://xenabrowser.net/?bookmark=98f9d901fdb2e95391fb4f5fdfac9097
```

### Como baixar pelo bookmark

1. Acesse o link acima (ou cole o código `98f9d901fdb2e95391fb4f5fdfac9097`
   no UCSC Xena em **File → Open bookmark**).
2. Selecione **File → Download**.
3. Salve como `liver.tsv` na raiz do projeto.

## Estrutura esperada do arquivo

O arquivo contém:

- **17 colunas de metadados** (`sample`, `sample_type`, `study`,
  `TCGA_GTEX_main_category`, `primary disease or tissue`, dados clínicos de
  sobrevida — `OS`, `DSS`, `DFI`, `PFI` — etc.);
- **212 colunas de expressão gênica** (genes da via KEGG `hsa05417`, escala
  log2).

| Grupo | Amostras |
|-------|----------|
| TCGA-LIHC (Primary Tumor) | 369 |
| GTEx Liver (Normal Tissue) | 110 |
| TCGA Solid Tissue Normal | 50 |
| TCGA Recurrent Tumor | 2 |
| **Total** | **531** |

> **Observação:** o pipeline atual compara **TCGA-LIHC (369)** × **GTEx Liver
> (110)**. As 50 amostras **TCGA Solid Tissue Normal** (tecido hepático normal
> adjacente ao tumor) e as 2 **Recurrent Tumor** estão disponíveis no arquivo e
> podem ser usadas em análises futuras (ex.: comparação pareada tumor × normal
> adjacente, sem efeito de lote TCGA/GTEx).
