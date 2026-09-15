# Nanotec 2026 — MEC e MMP9 no Hepatocarcinoma (HCC)

Análise transcriptômica (TCGA-LIHC + GTEx Liver, hub Toil via UCSCXenaTools) da
**remodelagem da matriz extracelular (MEC)** e da **superexpressão de MMP9** no
hepatocarcinoma, com foco em alvos do microambiente tumoral e perspectivas
translacionais para a nanobiotecnologia.

## Estrutura

```
nanotec2026/
├── dados/         (raw = download; processed = painel/transcriptoma filtrados)
├── scripts/       (00 download -> 01 auditoria -> 02 DE -> 03 sobrevivência -> 04 GSVA)
├── outputs/       (volcano, heatmap, tables, qc, survival, enrichment)
└── docs/          (RELATORIO_EXECUCAO.txt, DECISOES_TOMADAS.md, auditoria, manifesto)
```

## Como rodar (reproduzível)

```bash
cd nanotec2026
Rscript scripts/run_all.R     # roda todas as etapas na ordem
# ou, etapa a etapa:
#   Rscript scripts/00_download_data.R
#   Rscript scripts/01_audit_data.R
#   Rscript scripts/02_de_analysis.R
#   Rscript scripts/03_survival.R
#   Rscript scripts/04_enrichment.R
```

O download é feito com **UCSCXenaTools** (hub Toil). Os arquivos grandes
(transcriptoma completo) são ignorados pelo git e re-baixados automaticamente.

## Resultados principais

- **41/41 genes** do painel ECM/MMP/TIMP presentes; auditoria **aprovada**.
- **MMP9 up no HCC**: log2FC = +2.19, FDR = 8.5e-18 (limma, LIHC vs Normal).
- Também up: MMP1, MMP11, MMP12, MMP14, COL4A1/4A2, PLAU, FN1, VIM, SPP1.
- **MMP2 não diferencial** em bulk (log2FC +0.22) — resultado honesto documentado.
- Sobrevivência: MMP9 (HR 1.09, p=0.021), MMP1 (HR 1.17, p=6e-06) associados a
  pior sobrevida global.
- GSVA/ssGSEA: ECM-receptor interaction, EMT, Angiogenesis, Focal adhesion e
  Apical junction aumentados no tumor.

Leia a documentação completa em `docs/RELATORIO_EXECUCAO.txt` e
`docs/DECISOES_TOMADAS.md`.
