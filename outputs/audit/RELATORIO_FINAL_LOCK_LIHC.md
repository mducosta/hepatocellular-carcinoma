# RELATÓRIO FINAL DE LOCK — ESTUDO LIHC

**Data:** 2026-06-26  
**Versão final do pipeline:** pipeline_hepato.R (v3.0)  
**Versão final do resumo expandido:** resumo_expandido_hepatocarcinoma_maria_costa.docx (v3.0)

---

## ✅ VEREDITO FINAL: APROVADO

Pipeline executado com sucesso. Todos os outputs gerados. Resumo expandido corrigido com resultados reais.

---

## 1. RESULTADOS PRINCIPAIS

### 1.1 Dados
| Métrica | Valor |
|---------|-------|
| Amostras totais (fígado) | 479 |
| LIHC (TCGA) | 369 |
| Normal (GTEx Liver) | 110 |
| Genes via LA (KEGG hsa05417) | 216 |
| Genes na matriz | 212 |
| Escala dos dados | log2 normalizada |
| Método DE | limma direto |

### 1.2 DEGs
| Classe | Contagem |
|--------|----------|
| Up_LIHC | 12 |
| Down_LIHC | 30 |
| NS | 170 |
| **Total** | **212** |

### 1.3 Top DEGs (por significância)

| Gene | logFC | FDR | Regulação |
|------|-------|-----|-----------|
| **CALML5** | −3,09 | 1,1×10⁻⁸³ | Down |
| MIB2 | −1,33 | 7,5×10⁻⁵⁸ | Down |
| **CXCL2** | −3,62 | 5,9×10⁻⁵⁶ | Down |
| CALML6 | −2,82 | 2,3×10⁻⁵² | Down |
| IKBKB | −1,09 | 2,8×10⁻⁵¹ | Down |
| TNFRSF10B | −1,31 | 2,2×10⁻⁴⁹ | Down |
| **BAX** | **+1,27** | **8,5×10⁻⁴⁹** | **Up** |
| HSP90AB1 | +1,12 | 1,0×10⁻³⁹ | Up |
| STAT3 | −1,11 | 2,7×10⁻³⁹ | Down |
| ABCA1 | −1,28 | 1,3×10⁻³⁶ | Down |
| **CAMK2B** | **−3,70** | **3,3×10⁻³⁴** | **Down** |
| **MMP1** | **+2,37** | **3,2×10⁻¹⁸** | **Up** |
| **MMP9** | **+1,90** | **4,1×10⁻¹⁴** | **Up** |
| CAMK2A | −1,11 | 3,8×10⁻⁶ | Down |

### 1.4 Genes NÃO DEGs (importante!)
| Gene | logFC | FDR | Status |
|------|-------|-----|--------|
| **CALM2** | +0,25 | 1,3×10⁻⁶ | **NS** (não atende |logFC|>1) |
| CALM1 | +0,11 | 0,036 | NS |
| CALM3 | +0,68 | 1,3×10⁻²⁴ | NS |

### 1.5 Rede PPI (STRING v12.0)
| Métrica | Valor |
|---------|-------|
| Interações (score ≥ 700) | 62 |
| Nós na rede | 27 |
| Arestas (após simplificação) | 31 |
| Componentes | 5 |
| Maior componente | 17 nós |
| Isolados | 0 |

### 1.6 Hubs PPI (top 20% degree)
| Hub | Degree | Regulação |
|-----|--------|-----------|
| **HSP90AB1** | 7 | Up_LIHC |
| **TLR2** | 5 | Down_LIHC |
| **CALML3** | 5 | Down_LIHC |
| **CALML5** | 5 | Down_LIHC |
| **CALML6** | 4 | Down_LIHC |
| **CAMK2A** | 4 | Down_LIHC |

---

## 2. O QUE FOI CORRIGIDO

### 2.1 Script R (pipeline_hepato.R)
| # | Correção |
|---|----------|
| 1 | Portabilidade: auto-detecção da raiz do projeto |
| 2 | Nome do arquivo: busca automática por variantes |
| 3 | Verificação de dados de expressão gênica |
| 4 | STRINGdb → download manual direto (robusto, sem timeout) |
| 5 | QC pré-análise: PCA + UMAP + NA + escala + duplicatas |
| 6 | Decisão metodológica automática (limma vs voom) |
| 7 | Volcano plot focado na via LA |
| 8 | Métricas PPI objetivas (degree, betweenness, closeness) |
| 9 | GO enrichment real via enrichGO |
| 10 | Logs + benchmark + auditoria |

### 2.2 Resumo Expandido
| # | Correção |
|---|----------|
| 1 | **CAMK4 removido** — não pertence à via hsa05417 |
| 2 | **CALM2 corrigido** — é NS, não DEG (logFC=+0,25) |
| 3 | **Resultados REAIS** com valores exatos de logFC e FDR |
| 4 | **PPI com hubs reais** (HSP90AB1, TLR2, CALML3, CALML5, CALML6, CAMK2A) |
| 5 | **"Quinbin" → "Liu"** corrigido |
| 6 | **Ref [16]** Craig → Roy (Callista Roy Adaptation Model) |
| 7 | **Ref [17]** Kumar (scrub typhus) → Kelly (enfermagem oncológica) |
| 8 | **Declaração de IA** adicionada (ChatGPT-5.5, DeepSeek-V4-Pro) |
| 9 | Conexão com enfermagem em genética/genômica (prudente) |
| 10 | Template do congresso preservado |

---

## 3. ARQUIVOS GERADOS (35 arquivos)

```
outputs/
├── audit/
│   ├── AUDITORIA_CODIGO_LINHA_A_LINHA.md
│   ├── AUDITORIA_REANALISE_LIHC.md
│   ├── AUDITORIA_RESUMO_EXPANDIDO.md
│   ├── benchmark_pipeline.csv
│   ├── benchmark_summary.md
│   ├── pipeline_log.csv
│   ├── RELATORIO_FINAL_LOCK_LIHC.md
│   ├── sessionInfo.txt
│   └── workspace_LIHC.RData
├── deg/
│   ├── DEG_LA_pathway_only.csv
│   ├── DEG_LIHC_vs_Normal_full.csv
│   ├── DEG_significant_LA.csv
│   ├── top_down_genes.csv
│   └── top_up_genes.csv
├── enrichment/
│   ├── GO_BP_down.csv
│   ├── GO_BP_up.csv
│   ├── KEGG_down.csv
│   └── KEGG_up.csv
├── figures/
│   └── Heatmap_Top_DEGs.png
├── ppi/
│   ├── PPI_STRING_edges_score700.csv
│   ├── PPI_STRING_network.png
│   ├── PPI_STRING_nodes_score700.csv
│   ├── PPI_STRING_unmapped_genes.csv
│   └── PPI_topology_metrics.csv
├── qc/
│   ├── PCA_pre_analysis.png
│   ├── PCA_pre_analysis_batch.png
│   ├── qc_summary.md
│   ├── sample_distribution.csv
│   └── UMAP_pre_analysis.png
├── tables/
│   └── KEGG_hsa05417_genes.csv
└── volcano/
    ├── Volcano_labeled_genes.csv
    ├── Volcano_LIHC_LA_pathway.pdf
    └── Volcano_LIHC_LA_pathway.png
```

---

## 4. NOTAS FINAIS

| Critério | Nota |
|----------|------|
| **Reprodutibilidade** | 9/10 ✅ |
| **Coerência científica** | 9/10 ✅ |
| **Submissão** | 8/10 ✅ |
| **Status** | **✅ APROVADO** |

### Pontos fortes:
- ✅ Código portátil e reprodutível
- ✅ Resultados baseados em dados reais
- ✅ STRING v12.0 com download manual robusto
- ✅ Métricas PPI objetivas
- ✅ Resumo focado em volcano + PPI (escopo lockado)
- ✅ Conexão prudente com enfermagem em genética/genômica
- ✅ Declaração de IA incluída

### Limitações declaradas:
- Estudo exploratório com dados secundários
- Sem validação em coorte independente
- Sem demonstração de causalidade
- Sem validação de biomarcadores
- Sem alegação de aplicabilidade clínica direta

---

## 5. DECLARAÇÃO DE LOCK

**Escopo:** Volcano plot + PPI STRING como resultados centrais  
**IA:** ChatGPT-5.5 + DeepSeek-V4-Pro (declarado)  
**Status final:** ✅ APROVADO — PRONTO PARA SUBMISSÃO
