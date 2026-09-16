# Nanotec 2026 — MEC e MMP9 no Hepatocarcinoma (HCC)

Análise transcriptômica (TCGA-LIHC + GTEx Liver, hub Toil via UCSCXenaTools) da
**remodelagem da matriz extracelular (MEC)** e da **superexpressão de MMP9** no
hepatocarcinoma, com foco em alvos do microambiente tumoral e perspectivas
translacionais para a nanobiotecnologia.

## Estrutura

```
nanotec2026/
├── dados/         (raw = download; processed = painel/transcriptoma filtrados)
├── scripts/       (00 download -> 01 auditoria -> 02 DE -> 03 sobrevivência -> 04 GSVA -> 05 robustas)
├── outputs/       (volcano, heatmap, tables, qc, survival, enrichment, robust)
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
#   Rscript scripts/05_analises_robustas.R
```

O download é feito com **UCSCXenaTools** (hub Toil). Os arquivos grandes
(transcriptoma completo) são ignorados pelo git e re-baixados automaticamente.

## Resumo para submissão (Nanotec 2026 — Área temática 1: Biotecnologia e Nanobiotecnologia)

> **Título:** Remodelamento da matriz extracelular e superexpressão de MMP9 no
> hepatocarcinoma: assinatura do microambiente tumoral com valor prognóstico e
> perspectivas para a nanobiotecnologia
>
> **Autores:** Maria Eduarda Costa, Victória Oliveira Nascimento, Ryan de Paulo
> Santos, Heloisa Alves Guimarães
>
> **Resumo (193 palavras):** O hepatocarcinoma (HCC) progride em um microambiente
> tumoral (TME) fibrótico e imunossupressor, cuja matriz extracelular (MEC) densa
> constitui barreira à penetração de nanomateriais e reservatório de alvos
> moleculares. Este estudo objetivou caracterizar transcriptomicamente a
> reprogramação da MEC no HCC a partir de 529 amostras hepáticas (369 tumores
> TCGA-LIHC, 110 normais GTEx, 50 adjacentes). A expressão diferencial
> genoma-wide (limma) revelou 9.031 genes alterados, e o GSEA confirmou
> enriquecimento de transição epitelial-mesenquimal, resposta a interferons e
> apresentação de antígenos no tumor. No painel de 41 genes ECM/MMP/TIMP, MMP9
> (log2FC=+2,19; FDR=8,5e-18), MMP1, MMP11, MMP12, MMP14, COL4A1/4A2, PLAU, FN1,
> VIM e SPP1 mostraram-se superexpressos. A deconvolução celular (ssGSEA, atlas
> hepático Aizarani) evidenciou expansão de Treg, macrófagos M1, monócitos e
> endotélio, com MMP9 correlacionando-se a células de Kupffer (rho=0,54),
> fibroblastos/CAF e células estreladas. MMP9 covariou positivamente com
> checkpoints imunes (CTLA4 rho=0,69; PDCD1 rho=0,61; LAG3; PD-L1). Um escore
> prognóstico ECM-MMP (MMP9+MMP1+MMP12+MMP14) estratificou a sobrevida global
> (HR=1,62; p=3,5e-05), com efeito prognóstico de MMP9 predominante em homens
> (HR=1,13; p=0,01). Conclui-se que o eixo MMP9/MMPs-TIMPs define assinatura de
> TME imunologicamente ativo e fibrótico, apontando biomarcadores e alvos para
> nanosistemas responsivos a proteases e para superação de barreiras físicas à
> nanoterapia no HCC.

## Resultados principais

- **41/41 genes** do painel ECM/MMP/TIMP presentes; auditoria **aprovada**.
- **MMP9 up no HCC**: log2FC = +2.19, FDR = 8.5e-18 (limma, LIHC vs Normal).
- Também up: MMP1, MMP11, MMP12, MMP14, COL4A1/4A2, PLAU, FN1, VIM, SPP1.
- **MMP2 não diferencial** em bulk (log2FC +0.22) — resultado honesto documentado.
- Sobrevivência: MMP9 (HR 1.09, p=0.021), MMP1 (HR 1.17, p=6e-06) associados a
  pior sobrevida global.
- GSVA/ssGSEA: ECM-receptor interaction, EMT, Angiogenesis, Focal adhesion e
  Apical junction aumentados no tumor.

### Análises robustas complementares (Etapa 05 — `outputs/robust/`)

- **DE genoma-wide**: 9.031 DEGs (4.785 up / 4.246 down).
- **GSEA** (Hallmark/KEGG/GO:BP): EMT, resposta a interferon-α/γ, apresentação de
  antígenos e assinaturas proliferativas enriquecidas no tumor.
- **Deconvolução celular** (ssGSEA, atlas Aizarani): Treg, macrófagos M1,
  monócitos e endotélio aumentados; MMP9 correlaciona com Kupffer (ρ=0,54),
  fibroblastos/CAF (ρ=0,29) e células estreladas (ρ=0,22).
- **Checkpoints imunes**: MMP9 covariou com CTLA4 (ρ=0,69), TIM-3 (ρ=0,67),
  TIGIT (ρ=0,63), PD-1 (ρ=0,61), LAG3 e PD-L1.
- **Escore prognóstico ECM-MMP** (MMP9+MMP1+MMP12+MMP14): HR=1,62 (p=3,5e-05),
  log-rank p=0,005.
- **Sobrevida multi-desfecho** (OS/DSS/PFI/DFI) + Cox multivariado (ajustado por
  sexo): MMP9 OS HR=1,09 (p=0,021) estável no modelo ajustado.
- **Sexo**: MMP9 prognóstico em homens (HR=1,13, p=0,010), não em mulheres.

Leia a documentação completa em `docs/RELATORIO_EXECUCAO.txt` e
`docs/DECISOES_TOMADAS.md`.
