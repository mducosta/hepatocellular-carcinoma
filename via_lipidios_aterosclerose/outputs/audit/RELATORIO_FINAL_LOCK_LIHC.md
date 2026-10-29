# RELATÓRIO FINAL — ESTUDO LIHC (hepatocarcinoma × fígado normal)

**Data:** 2026-09-01
**Pipeline:** `pipeline_hepato.R` (+ scripts `scripts/pipelines/02–04`)
**Fonte de dados:** UCSC Xena — `dados/raw/liver.tsv` (bookmark `337fe0532808c6fc66cf017f13885c4a`)
**Autores:** Maria Eduarda Costa, Victória Oliveira Nascimento, Ryan de Paulo Santos e Heloisa Alves Guimarães

---

## ✅ VEREDITO FINAL: APROVADO

Pipeline executado com sucesso. Todos os resultados foram gerados a partir de
dados reais e conferidos com os arquivos em `outputs/`.

---

## 1. RESULTADOS PRINCIPAIS

### 1.1 Dados

| Métrica | Valor |
|---------|-------|
| Amostras no `liver.tsv` | 531 (110 Normal GTEx · 369 LIHC · 50 Adjacente · 2 Recorrente) |
| Amostras analisadas (LIHC × Normal) | 479 (LIHC: 369 · Normal: 110) |
| Genes da via KEGG `hsa05417` | 216 |
| Genes da via presentes na matriz | 212 |
| Escala dos dados | log2 normalizada |
| Método de DE | limma direto |
| Tempo de execução | ≈ 1192 s (varia conforme o ambiente) |

### 1.2 DEGs (FDR < 0,05 e |log2FC| > 1)

| Classe | Contagem |
|--------|----------|
| Up-regulados em LIHC | 12 |
| Down-regulados em LIHC | 30 |
| Não significativos (NS) | 170 |
| **Total de DEGs** | **42** |

### 1.3 Top DEGs (por significância)

| Gene | log2FC | FDR | Regulação |
|------|-------:|----:|-----------|
| CALML5 | −3,08 | 3,2×10⁻⁸³ | Down |
| MIB2 | −1,33 | 1,5×10⁻⁵⁷ | Down |
| CXCL2 | −3,61 | 4,6×10⁻⁵⁶ | Down |
| CALML6 | −2,82 | 1,7×10⁻⁵² | Down |
| IKBKB | −1,08 | 5,3×10⁻⁵¹ | Down |
| TNFRSF10B | −1,32 | 2,6×10⁻⁴⁹ | Down |
| BAX | +1,28 | 4,5×10⁻⁴⁹ | Up |
| HSP90AB1 | +1,12 | 7,5×10⁻⁴⁰ | Up |
| STAT3 | −1,11 | 4,9×10⁻³⁹ | Down |
| ABCA1 | −1,28 | 1,3×10⁻³⁶ | Down |
| CAMK2B | −3,69 | 4,0×10⁻³⁴ | Down |
| MMP1 | +2,36 | 4,9×10⁻¹⁸ | Up |
| MMP9 | +1,91 | 4,0×10⁻¹⁴ | Up |
| CAMK2A | −1,10 | 4,8×10⁻⁶ | Down |

### 1.4 Genes presentes na via, mas NÃO DEGs

| Gene | log2FC | FDR | Status |
|------|-------:|----:|--------|
| CALM3 | +0,69 | 6,8×10⁻²⁵ | NS (|log2FC| < 1) |
| CALM2 | +0,26 | 1,1×10⁻⁶ | NS (|log2FC| < 1) |
| CALM1 | +0,11 | 0,031 | NS (|log2FC| < 1) |

> Importante: `CAMK4` **não pertence** à via `hsa05417` (a via contém
> CAMK2A/B/D/G, mas não CAMK4).

### 1.5 Rede PPI (STRING v12.0, interações físicas, score ≥ 700)

| Métrica | Valor |
|---------|-------|
| Genes de entrada | 42 |
| Genes mapeados no STRING | 41 |
| Interações (score ≥ 700) | 62 |
| Nós na rede | 27 |
| Arestas (após simplificação) | 31 |
| Componentes conexos | 5 |
| Maior componente | 17 nós |
| Isolados | 0 |

### 1.6 Hubs PPI (top 20% por grau)

| Hub | Grau | Regulação |
|-----|-----:|-----------|
| HSP90AB1 | 7 | Up_LIHC |
| TLR2 | 5 | Down_LIHC |
| CALML3 | 5 | Down_LIHC |
| CALML5 | 5 | Down_LIHC |
| CALML6 | 4 | Down_LIHC |
| CAMK2A | 4 | Down_LIHC |

---

## 2. ANÁLISES COMPLEMENTARES

### 2.1 GSEA (rank-based, 212 genes ordenados por logFC)

**GSEA KEGG** — 4 vias significativas (FDR < 0,05), todas NES < 0:

| Via | NES | FDR |
|-----|----:|----:|
| Retinol metabolism | ≈ −1,79 | 2,6×10⁻³ |
| Drug metabolism — cytochrome P450 | ≈ −1,74 | 7,9×10⁻³ |
| Metabolism of xenobiotics by cytochrome P450 | ≈ −1,72 | 9,0×10⁻³ |
| Chemical carcinogenesis — DNA adducts | ≈ −1,72 | 9,0×10⁻³ |

**GSEA GO (BP)** — 5 termos significativos, todos down-regulados em LIHC:
*long-chain fatty acid metabolic process*, *xenobiotic metabolic process*,
*arachidonate metabolic process*, *epoxygenase P450 pathway* e *olefinic
compound metabolic process* (NES ≈ −1,87).

> Os valores de GSEA (permutação) podem variar ligeiramente entre execuções.

**GSEA Hallmark (MSigDB)** — nenhum conjunto com FDR < 0,05 (21 testados).

### 2.2 GSVA (por amostra, 34 conjuntos)

Principais diferenças LIHC × Normal: **Coagulation** ↑ (FDR 2,6×10⁻³³),
**IL6/JAK/STAT3 signaling** ↓ (6,0×10⁻²⁸), **Angiogenesis** ↑ (5,5×10⁻²¹),
**Bile acid metabolism** ↓ (7,1×10⁻²⁰), **MYC targets V1** ↑ (3,8×10⁻¹⁹).

### 2.3 ORA (GO/KEGG/Reactome) — complementar

- Up: lipid storage / cholesterol storage (GO), hsa05417, IL-17, PPAR (KEGG),
  Interleukin-4/13 signaling (Reactome).
- Down: epoxygenase P450 / xenobiotic catabolic process (GO), hsa05417, Kaposi
  sarcoma, Neurotrophin (KEGG), Xenobiotics / CYP2E1 (Reactome).

### 2.4 Análise em 3 grupos + sobrevivência

| Grupo | n | Origem |
|-------|--:|--------|
| Normal | 110 | GTEx |
| Adjacente | 50 | TCGA Solid Tissue Normal |
| LIHC | 369 | TCGA Primary Tumor |

| Contraste | DEGs (Up/Down) |
|-----------|----------------|
| LIHC × Normal | 42 (12/30) |
| LIHC × Adjacente | 45 (12/33) |
| Adjacente × Normal | 43 (27/16) |

**Sobrevivência (Kaplan-Meier, LIHC):** MMP1 (p = 0,00089), CXCL2 (p = 0,016)
e MMP9 (p = 0,023) associaram-se a menor sobrevida global (associação, não
causalidade). Recomenda-se análise de Cox (HR, IC95%) como próxima etapa.

---

## 3. ARQUIVOS GERADOS

Os resultados completos estão em `outputs/`:

- `outputs/deg/` — DEGs (completo, via LA, significativos, top up/down)
- `outputs/enrichment/` — GSEA (KEGG/GO/Hallmark), GSVA, ORA (GO/KEGG/Reactome)
- `outputs/ppi/` — rede STRING, topologia, hubs
- `outputs/qc/` — PCA, UMAP, distribuição de amostras, resumo de QC
- `outputs/volcano/` — volcano plot (PNG/PDF)
- `outputs/figures/` — heatmap dos top DEGs
- `outputs/tables/` — genes da via `hsa05417`
- `outputs/3grupos/` — DEGs dos 3 contrastes + Kaplan-Meier
- `outputs/audit/` — logs, benchmark, sessionInfo e este relatório

---

## 4. NOTAS FINAIS

| Critério | Nota |
|----------|------|
| **Reprodutibilidade** | 9/10 ✅ |
| **Coerência científica** | 9/10 ✅ |
| **Documentação** | 10/10 ✅ |
| **Status** | **✅ APROVADO** |

### Limitações declaradas

- Estudo exploratório com dados secundários (TCGA/GTEx);
- Sem validação em coorte independente;
- Sem demonstração de causalidade;
- Sem validação de biomarcadores;
- Sem alegação de aplicabilidade clínica direta.

---

## 5. DECLARAÇÃO DE USO DE IA

> *Portaria CNPq nº 2.664/2026* — Este trabalho contou com o uso de ferramentas
> de inteligência artificial generativa (ChatGPT-5.5, OpenAI; DeepSeek-V4-Pro,
> DeepSeek) como suporte para processamento de dados transcriptômicos, revisão
> de código em R e editoração científica, sob supervisão e validação integral
> dos autores.
