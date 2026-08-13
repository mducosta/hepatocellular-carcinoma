# AUDITORIA DE REANÁLISE LIHC — RESULTADOS FINAIS

**Data:** 2026-06-26  
**Status:** ✅ APROVADO  
**Pipeline:** pipeline_hepato.R (v3.0, portátil, reprodutível)

---

## 1. EXECUÇÃO

| Etapa | Status | Detalhe |
|-------|--------|---------|
| Detecção do projeto | ✅ | C:/Users/oorie/OneDrive/Documentos/TRABALHOS/LIVER |
| Arquivo de dados | ✅ | liver_lip_aterosclerose.tsv (531 × 230) |
| Metadados | ✅ | 17 colunas clínicas |
| Expressão gênica | ✅ | 212 genes da via hsa05417 |
| Amostras hepáticas | ✅ | 479 (369 LIHC, 110 Normal) |
| QC (PCA, UMAP) | ✅ | Gerados |
| Genes KEGG | ✅ | 216 genes (212 na matriz) |
| Análise diferencial | ✅ | limma direto (log2 normalizada) |
| DEGs | ✅ | 42 (12 Up, 30 Down) |
| Volcano plot | ✅ | PNG + PDF gerados |
| PPI STRING | ✅ | 62 interações, 6 hubs |
| Enriquecimento KEGG (ORA) | ✅ | 47 vias (Up), 107 vias (Down) |
| GO BP (ORA) | ✅ | 78 termos (Up), 44 termos (Down) |
| Heatmap | ✅ | Top 42 DEGs |
| Benchmark | ✅ | 676.92 segundos |

---

## 2. VERIFICAÇÕES DE QUALIDADE

| Verificação | Resultado |
|-------------|-----------|
| Valores NA | 0 (0,00%) ✅ |
| Escala dos dados | log2 [0,00 – 21,37] ✅ |
| Integer-like | FALSE ✅ |
| Genes duplicados | 0 ✅ |
| Amostras duplicadas | 0 ✅ |
| Genes zero-var | 0 ✅ |
| Batch TCGA/GTEx | Identificado (PCA mostra separação por estudo) ⚠️ |
| PCA PC1 | 12,4% |
| PCA PC2 | 11,5% |
| Decisão metodológica | limma direto (log2) ✅ |

---

## 3. GENES DA VIA vs RESULTADOS

### 3.1 Genes citados no resumo × status real

| Gene | Na via? | DEG? | logFC | FDR | Direção |
|------|---------|------|-------|-----|---------|
| CALML5 | ✅ | ✅ | −3,09 | 1,1e-83 | Down |
| CXCL2 | ✅ | ✅ | −3,62 | 5,9e-56 | Down |
| BAX | ✅ | ✅ | +1,27 | 8,5e-49 | Up |
| CAMK2B | ✅ | ✅ | −3,70 | 3,3e-34 | Down |
| MMP1 | ✅ | ✅ | +2,37 | 3,2e-18 | Up |
| MMP9 | ✅ | ✅ | +1,90 | 4,1e-14 | Up |
| CAMK2A | ✅ | ✅ | −1,11 | 3,8e-06 | Down |
| **CALM2** | ✅ | **❌ NS** | +0,25 | 1,3e-06 | — |
| **CALM1** | ✅ | **❌ NS** | +0,11 | 0,036 | — |
| **CALM3** | ✅ | **❌ NS** | +0,68 | 1,3e-24 | — |
| **CAMK4** | **❌ NÃO** | — | — | — | — |

### 3.2 Conclusões da verificação
- **CALM2**: presente na via mas NÃO é DEG (|logFC| < 1). O resumo original erroneamente o tratava como DEG principal.
- **CAMK4**: NÃO pertence à via hsa05417. O resumo original o citava incorretamente.
- **CALML5**: confirmado como gene mais significativo (FDR = 1,1×10⁻⁸³), corroborando parcialmente o resumo original.
- **CAMK2B**: correto mencionar (pertence à via e é DEG), mas o resumo original citava CAMK4 em seu lugar.

---

## 4. REDE PPI — HUBS CONFIRMADOS

| Hub | Degree | Betweenness | Função |
|-----|--------|-------------|--------|
| HSP90AB1 | 7 | 0,2800 | Chaperona, estabilidade proteica |
| TLR2 | 5 | 0,1262 | Receptor imune inato |
| CALML3 | 5 | 0,0200 | Calmodulina-like |
| CALML5 | 5 | 0,0200 | Calmodulina-like |
| CALML6 | 4 | 0,0369 | Calmodulina-like |
| CAMK2A | 4 | 0,0000 | Quinase dependente de Ca²⁺/calmodulina |

---

## 5. CORREÇÕES APLICADAS AO RESUMO

1. ❌ CAMK4 → ✅ CAMK2A/CAMK2B
2. ❌ CALM2 como DEG → ✅ CALM2 é NS
3. ❌ Dados inventados → ✅ Resultados reais com valores exatos
4. ❌ "Quinbin" → ✅ "Liu"
5. ❌ Ref [16] Craig → ✅ Roy (Callista Roy)
6. ❌ Ref [17] scrub typhus → ✅ Kelly (enfermagem oncológica)
7. ❌ Alegações prognósticas → ✅ Linguagem exploratória prudente
8. ➕ Declaração de IA adicionada

---

## 6. CONFORMIDADE COM O LOCK

| Requisito | Status |
|-----------|--------|
| Volcano plot como resultado central | ✅ |
| PPI STRING real como resultado central | ✅ |
| PCA/UMAP apenas como QC | ✅ |
| Enriquecimento como complementar | ✅ |
| Sem alegação de biomarcador | ✅ |
| Sem alegação de alvo terapêutico | ✅ |
| Sem alegação de prognóstico | ✅ |
| Sem alegação de decisão clínica direta | ✅ |
| Conexão com enfermagem genética/genômica | ✅ |
| Template do congresso preservado | ✅ |
| Declaração de IA | ✅ |
