# Vias para análise futura (KEGG + Reactome)

Lista de vias candidatas para ampliar a análise de enriquecimento funcional
(GSEA/ORA) do estudo LIHC × Normal. As vias estão separadas por prioridade.

---

## ⭐ VIA PRINCIPAL (próxima a ser testada)

| Via | ID | Justificativa |
|-----|----|---------------|
| **Drug metabolism — cytochrome P450** | **KEGG `hsa00982`** | **Via principal.** Foi a 2ª via mais significativa no GSEA KEGG (NES = −1,73; FDR = 5,5×10⁻³), com core enrichment dos genes CYP2C8, CYP2C9, CYP2B6, CYP2A6, CYP2A7 (down-regulados em LIHC). Coerente com a ORA (epoxygenase P450, xenobióticos) e com a Reactome (Xenobiotics, CYP2E1 reactions). É o eixo metabólico mais forte do estudo. |

**Objetivo do teste futuro:** caracterizar em detalhe o eixo de metabolismo de
xenobióticos/CYP450 no hepatocarcinoma — verificar se a perda de expressão de
CYP450 correlaciona com prognóstico, etiologia (NASH/álcool/hepatite) e resposta
a fármacos (sorafenibe etc.).

---

## KEGG — vias candidatas

| Via | ID | Justificativa |
|-----|----|---------------|
| **Drug metabolism — cytochrome P450** ⭐ | hsa00982 | Via principal (ver acima) |
| Hepatocellular carcinoma | hsa05225 | Via específica do CHC — obrigatória |
| Metabolism of xenobiotics by cytochrome P450 | hsa00980 | Coerente com hsa00982 (GSEA significativa, NES −1,71) |
| Retinol metabolism | hsa00830 | GSEA significativa (NES −1,84) — metabolismo lipídico |
| Chemical carcinogenesis — DNA adducts | hsa05204 | GSEA significativa (NES −1,71) |
| NAFLD | hsa04932 | Contexto etiológico (já enriquecida na ORA down) |
| Alcoholic liver disease | hsa04936 | Contexto etiológico (já enriquecida na ORA down) |
| Complement and coagulation cascades | hsa04610 | Coerente com GSVA (coagulação ↑, FDR 2,6×10⁻³³) |
| PPAR signaling pathway | hsa03320 | Enriquecida nos genes up |
| PI3K-Akt signaling pathway | hsa04151 | Enriquecida — proliferação |
| MAPK signaling pathway | hsa04010 | Enriquecida nos genes down |
| IL-17 signaling pathway | hsa04657 | Enriquecida nos genes up |
| Lipid and atherosclerosis | hsa05417 | Via original do estudo (referência) |

---

## Reactome — vias candidatas

| Pathway | ID | Justificativa |
|---------|----|---------------|
| Metabolism of lipids | R-HSA-556833 | Metabolismo lipídico (contexto da via) |
| Fatty acid metabolism | R-HSA-8978868 | Ácidos graxos (CYP/ACOX) |
| Cytochrome P450 — arranged by substrate type | R-HSA-211897 | **Já enriquecida (ORA down, FDR 6,8×10⁻⁶)** |
| Xenobiotics | R-HSA-211981 | **Já enriquecida (ORA down, FDR 2,5×10⁻⁸)** |
| Signaling by Interleukins | R-HSA-449147 | **Já enriquecida (ORA up)** — inflamação |
| Toll-like Receptor Cascades | R-HSA-168898 | **Já enriquecida (ORA down)** — imunidade inata |
| Extracellular matrix organization | R-HSA-1474244 | MMP1/MMP9 — invasão |
| Regulation of lipid metabolism by PPARα | R-HSA-400206 | PPAR — metabolismo lipídico |
| NLRP3 inflammasome | R-HSA-844456 | **Já enriquecida (ORA up)** — piroptose |

---

## Observações

- As vias marcadas com "já enriquecida" foram testadas pelo script
  `scripts/analise_reactome.R` (resultados em `results/enrichment/Reactome_*`).
- Para aprofundar a **via principal `hsa00982`**, sugere-se: baixar o
  **transcriptoma completo** (ver `docs/tutorial_dados_xena.md`) para que o
  GSEA/ORA tenha poder estatístico adequado.
