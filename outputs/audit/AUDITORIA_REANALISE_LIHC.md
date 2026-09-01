# AUDITORIA DE REANÁLISE LIHC — RESULTADOS FINAIS

**Data:** 2026-09-01
**Status:** ✅ APROVADO
**Pipeline:** `pipeline_hepato.R` (+ scripts `scripts/pipelines/02–04`)
**Autores:** Maria Eduarda Costa, Victória Oliveira Nascimento, Ryan de Paulo Santos e Heloisa Alves Guimarães

---

## 1. EXECUÇÃO

| Etapa | Status | Detalhe |
|-------|--------|---------|
| Detecção do projeto | ✅ | Automática (raiz do repositório) |
| Arquivo de dados | ✅ | `dados/raw/liver.tsv` (531 amostras × 212 genes da via) |
| Metadados | ✅ | 17 colunas clínicas/identificadoras |
| Expressão gênica | ✅ | 212 genes da via KEGG `hsa05417` |
| Amostras hepáticas (LIHC × Normal) | ✅ | 479 (369 LIHC, 110 Normal) |
| QC (PCA, UMAP) | ✅ | Gerados |
| Genes KEGG | ✅ | 216 genes (212 na matriz) |
| Análise diferencial | ✅ | limma direto (log2 normalizada) |
| DEGs | ✅ | 42 (12 Up, 30 Down) |
| Volcano plot | ✅ | PNG + PDF gerados |
| PPI STRING | ✅ | 62 interações, 6 hubs |
| Enriquecimento KEGG (ORA) | ✅ | 47 vias (Up), 107 vias (Down) |
| GO BP (ORA) | ✅ | 78 termos (Up), 44 termos (Down) |
| Reactome (ORA) | ✅ | 13 termos (Up), 32 termos (Down) |
| GSEA (KEGG/GO/Hallmark) | ✅ | 4 / 5 / 0 vias significativas |
| GSVA | ✅ | 34 conjuntos |
| Heatmap | ✅ | Top 42 DEGs |
| Sobrevivência (Kaplan-Meier) | ✅ | 9 genes testados |

---

## 2. VERIFICAÇÕES DE QUALIDADE

| Verificação | Resultado |
|-------------|-----------|
| Valores NA | 0 (0,00%) ✅ |
| Escala dos dados | log2 [0,00 – 21,37] ✅ |
| Integer-like | FALSE ✅ |
| Genes duplicados | 0 ✅ |
| Amostras duplicadas | 0 ✅ |
| Genes zero-variância | 0 ✅ |
| Batch TCGA/GTEx | Identificado (PCA mostra separação por estudo) ⚠️ |
| Decisão metodológica | limma direto (log2) ✅ |

---

## 3. GENES DA VIA × STATUS REAL

### 3.1 Genes citados × status real

| Gene | Na via? | DEG? | log2FC | FDR | Direção |
|------|---------|------|-------:|----:|---------|
| CALML5 | ✅ | ✅ | −3,08 | 3,2×10⁻⁸³ | Down |
| CXCL2 | ✅ | ✅ | −3,61 | 4,6×10⁻⁵⁶ | Down |
| BAX | ✅ | ✅ | +1,28 | 4,5×10⁻⁴⁹ | Up |
| CAMK2B | ✅ | ✅ | −3,69 | 4,0×10⁻³⁴ | Down |
| MMP1 | ✅ | ✅ | +2,36 | 4,9×10⁻¹⁸ | Up |
| MMP9 | ✅ | ✅ | +1,91 | 4,0×10⁻¹⁴ | Up |
| CAMK2A | ✅ | ✅ | −1,10 | 4,8×10⁻⁶ | Down |
| CALM3 | ✅ | ❌ NS | +0,69 | 6,8×10⁻²⁵ | — |
| CALM2 | ✅ | ❌ NS | +0,26 | 1,1×10⁻⁶ | — |
| CALM1 | ✅ | ❌ NS | +0,11 | 0,031 | — |
| **CAMK4** | **❌ NÃO** | — | — | — | — |

### 3.2 Conclusões da verificação

- **CALM2** está na via, mas **não é DEG** (|log2FC| < 1). Não deve ser tratado
  como DEG principal.
- **CAMK4** **não pertence** à via `hsa05417` (a via contém CAMK2A/B/D/G, não
  CAMK4). Não deve ser citado como resultado da via.
- **CALML5** confirmado como o gene mais significativo (FDR = 3,2×10⁻⁸³).

---

## 4. REDE PPI — HUBS CONFIRMADOS

| Hub | Grau | Função |
|-----|-----:|--------|
| HSP90AB1 | 7 | Chaperona, estabilidade proteica |
| TLR2 | 5 | Receptor imune inato |
| CALML3 | 5 | Calmodulina-like |
| CALML5 | 5 | Calmodulina-like |
| CALML6 | 4 | Calmodulina-like |
| CAMK2A | 4 | Quinase dependente de Ca²⁺/calmodulina |

---

## 5. CORREÇÕES APLICADAS AO RESUMO EXPANDIDO

1. ❌ CAMK4 → ✅ CAMK2A/CAMK2B (CAMK4 não pertence à via)
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
| Sem alegação de prognóstico (apenas associação) | ✅ |
| Sem alegação de decisão clínica direta | ✅ |
| Conexão com enfermagem genética/genômica | ✅ |
| Declaração de IA | ✅ |
