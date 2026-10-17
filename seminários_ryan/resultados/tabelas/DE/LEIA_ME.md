# Resultados de expressão diferencial

Cada comparação gera:

- `DEG_<comparacao>_completo.csv` — todos os genes testados.
- `DEG_<comparacao>_significativos.csv` — apenas FDR < 0,05 e |logFC| > 1.
- `TOP20_<comparacao>.csv` — top 20 up + top 20 down.

## Colunas

| Coluna | Descrição |
|---|---|
| `gene` | símbolo HUGO |
| `comparison` | nome da comparação |
| `direction` | direção do contraste |
| `logFC` | log2 fold change |
| `AveExpr` | expressão média (log2) |
| `t` | t moderada do limma |
| `P.Value` | p nominal |
| `adj.P.Val` | FDR (BH) |
| `B` | log-odds bayesiano |
| `regulation` | Up / Down / NS |

## Arquivos-resumo

- `RESUMO_DEGs_TODAS_COMPARACOES.csv` — contagens up/down/total por comparação.
- `RESUMO_TOP5_GENES.csv` — top 5 up e down de cada comparação.
- `COMPARACAO_ComBat_antes_depois.csv` — efeito da demonstração do ComBat.

## Comparações

| Nome | Contraste |
|---|---|
| `TCGA_Tumor_vs_TCGA_Adjacent` | tumor − adjacente (TCGA) |
| `TCGA_Tumor_vs_Adjacent_PAREADO` | tumor − adjacente, paciente como bloqueio |
| `TCGA_Tumor_vs_GTEx_Normal` | tumor TCGA − fígado GTEx |
| `TCGA_Adjacent_vs_GTEx_Normal` | adjacente TCGA − fígado GTEx |
| `GTEx_Male_vs_Female` | masculino − feminino (GTEx) |
| `GTEx_Age60plus_vs_lt60` | ≥60 − <60 (GTEx) |
| `TCGA_Stage_III_IV_vs_I_II` | III–IV − I–II (TCGA tumor) |
| `TCGA_Grade_G3_G4_vs_G1_G2` | G3–G4 − G1–G2 (TCGA tumor) |
| `TCGA_Male_vs_Female` | masculino − feminino (TCGA tumor) |
| `TCGA_Age65plus_vs_lt65` | ≥65 − <65 (TCGA tumor) |
| `LApathway_*` | via KEGG hsa05417 |
| `*_ComBat*` | demonstração de correção de lote |
