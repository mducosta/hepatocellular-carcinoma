# Hepatocarcinoma — Estudos Transcriptômicos

Repositório que reúne três estudos independentes sobre o transcriptoma do
**hepatocarcinoma (HCC)**, usando dados públicos de RNA-seq do **TCGA-LIHC** e
do **fígado normal do GTEx** (hub Toil, UCSC Xena), além de coortes de
validação do GEO.

Cada subpasta é um projeto autocontido, com dados processados, scripts,
resultados, documentação e auditorias próprios.

---

## 1. Subprojetos

| Pasta | Estudo | Abordagem principal | Entrada | Documentação |
|---|---|---|---|---|
| [`via_lipidios_aterosclerose/`](via_lipidios_aterosclerose/) | Via KEGG `hsa05417` (*Lipid and Atherosclerosis*) no LIHC | Expressão diferencial LIHC × fígado normal; QC, rede PPI/STRING, GSEA, GSVA, ORA e Reactome | `pipeline_hepato.R` | [`README.md`](via_lipidios_aterosclerose/README.md) · [`paper_LIHC.pdf`](via_lipidios_aterosclerose/paper_LIHC.pdf) |
| [`analise_covariaveis_deconvolucao/`](analise_covariaveis_deconvolucao/) | Transcriptoma do LIHC ajustado por clínica, lote e composição celular | Tumor × adjacente × GTEx; covariáveis clínicas, ComBat, deconvolução (EPIC, quanTIseq, xCell), sensibilidades e validação externa | `scripts/executar_todas_as_etapas.R` | [`LEIA_ME.md`](analise_covariaveis_deconvolucao/LEIA_ME.md) · [artigo final (DOCX)](analise_covariaveis_deconvolucao/artigo_final/artigo_expressao_diferencial_covariaveis_deconvolucao.docx) |
| [`nanotec2026/`](nanotec2026/) | Remodelagem da matriz extracelular e MMP9 no HCC | Assinatura do microambiente tumoral, sobrevivência, GSVA/ssGSEA e análises robustas | `scripts/run_all.R` | [`README.md`](nanotec2026/README.md) · [`docs/DECISOES_TOMADAS.md`](nanotec2026/docs/DECISOES_TOMADAS.md) |

### 1.1 Resumo dos achados principais

- **Via `hsa05417` (Lipid and Atherosclerosis):** 42 DEGs (12 up / 30 down) na
  comparação LIHC × normal, com perda de genes de metabolismo de xenobióticos
  (CYP2C8, CYP2C9, CYP2B6, CYP2A6/7) e ganho de genes inflamatórios e de
  remodelamento (MMP1, MMP9, CD36, VCAM1). Resultados detalhados em
  [`via_lipidios_aterosclerose/README.md`](via_lipidios_aterosclerose/README.md).
- **Covariáveis e deconvolução:** o eixo tumor × adjacente é o mais robusto
  (4.130 DEGs no basal; 3.710 no pareado; 3.740 no modelo clínico), enquanto a
  assinatura de estágio se torna nula após ajuste por centro/composição e a de
  grau se atenua fortemente. A comparação TCGA × GTEx é estruturalmente
  confundida com estudo/lote. Síntese em
  [`analise_covariaveis_deconvolucao/LEIA_ME.md`](analise_covariaveis_deconvolucao/LEIA_ME.md).
- **Nanotec 2026:** MMP9 está superexpresso no HCC (log2FC = +2,19;
  FDR = 8,5 × 10⁻¹⁸) e compõe, com MMP1/MMP12/MMP14, um escore prognóstico
  associado a pior sobrevida global (HR = 1,62; p = 3,5 × 10⁻⁵). Resumo completo
  em [`nanotec2026/README.md`](nanotec2026/README.md).

---

## 2. Estrutura do repositório

```text
hepatocellular-carcinoma/
├── README.md                              # Este índice
├── LICENSE                                # MIT
├── .gitignore                             # Dados grandes e artefatos locais
│
├── via_lipidios_aterosclerose/            # Estudo 1 — via KEGG hsa05417
│   ├── README.md
│   ├── requirements.txt
│   ├── pipeline_hepato.R
│   ├── paper_LIHC.pdf
│   ├── scripts/
│   ├── dados/
│   ├── outputs/
│   └── docs/
│
├── analise_covariaveis_deconvolucao/      # Estudo 2 — covariáveis, lote e composição
│   ├── LEIA_ME.md
│   ├── relatorio_consolidado_seminario.md
│   ├── scripts/
│   ├── dados/
│   ├── resultados/
│   ├── auditorias/
│   ├── docs/
│   └── artigo_final/
│
└── nanotec2026/                           # Estudo 3 — MEC e MMP9
    ├── README.md
    ├── scripts/
    ├── dados/
    ├── outputs/
    └── docs/
```

---

## 3. Dados e reprodutibilidade

- **Fontes públicas:** UCSC Xena/Toil (TCGA-LIHC e GTEx), GDC e GEO
  (GSE14520, GSE36376, GSE76427 e GSE62232). Nenhum dado identificável de
  paciente é usado.
- **Arquivos grandes não versionados:** transcriptomas brutos, contagens do GDC,
  matrizes de expressão das coortes externas e tabelas derivadas de grande porte
  são ignorados pelo Git e podem ser obtidos ou regenerados com os scripts de
  cada projeto.
- **Dados leves versionados:** matrizes processadas, metadados, tabelas-resumo,
  figuras e relatórios necessários para inspecionar os resultados.
- **Ambientes:** R ≥ 4.0 (os estudos avançados usam R 4.6/Bioconductor 3.23) e
  Python 3. As versões exatas estão em `sessionInfo` e arquivos de versões
  dentro de cada subprojeto.

### 3.1 Execução rápida

```bash
# Estudo 1 — via hsa05417
cd via_lipidios_aterosclerose
Rscript pipeline_hepato.R

# Estudo 2 — pipeline completo (covariáveis, ComBat, deconvolução, validações)
cd ../analise_covariaveis_deconvolucao
Rscript scripts/executar_todas_as_etapas.R

# Estudo 3 — MEC e MMP9
cd ../nanotec2026
Rscript scripts/run_all.R
```

Consulte o README/LEIA_ME de cada pasta antes de executar: alguns passos
dependem de downloads grandes e de pacotes específicos.

---

## 4. Autores

- **Maria Eduarda Costa**
- **Victória Oliveira Nascimento**
- **Ryan de Paulo Santos**
- **Heloisa Alves Guimarães**

---

## 5. Licença

[MIT](LICENSE) © 2026 Maria Eduarda Costa, Victória Oliveira Nascimento, Ryan
de Paulo Santos e Heloisa Alves Guimarães.

---

## 6. Declaração de uso de inteligência artificial

Em conformidade com a Portaria CNPq nº 2.664/2026, declara-se que ferramentas
de IA generativa foram utilizadas como assistentes de programação, organização
do pipeline, revisão de código e redação técnica. Todas as rotinas foram
executadas, conferidas e validadas contra os arquivos de saída; a interpretação
dos resultados, a definição do desenho analítico e a responsabilidade
intelectual pelo conteúdo permanecem sob responsabilidade dos autores. Nenhum
dado identificável de paciente foi usado.
