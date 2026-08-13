# Implicações clínicas para a Enfermagem

> Documento de apoio à interpretação clínica dos resultados transcriptômicos do
> estudo **LIHC (hepatocarcinoma) × fígado normal (GTEx)** — via *Lipid and
> Atherosclerosis* (KEGG `hsa05417`).

---

## 1. Resumo dos achados com relevância clínica

| Achado molecular | Evidência | Implicação clínica |
|------------------|-----------|--------------------|
| **Coagulação aumentada no tumor** | GSVA: HALLMARK_COAGULATION ↑ (FDR 2,6×10⁻³³) | Risco de tromboembolismo venoso (TEV) |
| **Angiogênese e proliferação ↑** | GSVA: Angiogenesis ↑, MYC targets ↑ | Agressividade tumoral; base de anti-angiogênicos |
| **Perda de metabolismo de xenobióticos/CYP450** | GSEA: hsa00982 (FDR 5,5×10⁻³); ORA: epoxygenase P450 | Perda de função hepatocitária; cuidado com fármacos metabolizados no fígado |
| **Inflamação (IL-4/IL-13, TLR) ↑** | Reactome ORA up: interleucinas; KEGG: IL-17 | Processo inflamatório crônico no tumor |
| **Remodelamento de matriz (MMP1/MMP9) ↑** | DEGs up; **associados a pior sobrevida** | Invasão/metástase; prognóstico |
| **Bile acid metabolism ↓** | GSVA: HALLMARK_BILE_ACID_METABOLISM ↓ | Disfunção hepatobiliar |

---

## 2. Como o enfermeiro pode aproveitar (por eixo)

### 🩸 2.1 Vigilância de tromboembolismo (coagulação ↑)
- Aplicar **escala de risco de TEV** (ex.: Khorana) na admissão do paciente com CHC.
- **Educar o paciente e a família** sobre sinais de TVP (dor/edema unilateral de membro) e TEP (dispneia súbita, dor torácica, taquicardia).
- Estimular **mobilização precoce** e medidas mecânicas (meias de compressão) quando indicado.
- Registrar e comunicar ao médico qualquer suspeita — TEV é uma das principais causas evitáveis de morte no paciente oncológico.

### 🔥 2.2 Manejo da inflamação e sintomas sistêmicos
- O tumor apresenta **assinatura inflamatória** (interleucinas, TLR) — base para compreender sintomas como **febre, fadiga, caquexia** e resposta inflamatória sistêmica.
- Avaliar sistematicamente **estado nutricional e perda de peso** (caquexia associada ao câncer).
- Correlacionar com o **controle de comorbidades inflamatórias** (hepatite, esteatose).

### 🫀 2.3 Educação sobre fatores de risco metabólicos
- Os resultados ligam o CHC a **metabolismo lipídico** (CD36 ↑, ABCA1 ↓, bile acids ↓) — reforça o aconselhamento sobre:
  - **Álcool** (fator de risco para CHC e doença hepática);
  - **Síndrome metabólica / esteatose hepática (NAFLD/NASH)**;
  - **Hepatites virais** (B e C) — adesão à vacinação e ao tratamento antiviral;
  - **Dieta e controle de lipídios/glicemia**.

### 💊 2.4 Segurança farmacológica (CYP450 ↓)
- A **perda de expressão de CYP450** no tumor sugere alteração do metabolismo hepático de fármacos.
- O enfermeiro deve **vigiar toxicidade medicamentosa** (ex.: anticoagulantes, analgésicos, antineoplásicos orais metabolizados no fígado) e relatar interações.
- Reforçar a **reconciliação medicamentosa** e a comunicação com a farmácia clínica.

### 📈 2.5 Prognóstico e educação do paciente
- **MMP1, CXCL2 e MMP9** (up-regulados) associaram-se a **pior sobrevida global** (Kaplan-Meier, p < 0,05).
- Esse tipo de assinatura molecular pode, no futuro, apoiar a **estratificação de risco** — cabe ao enfermeiro compreender e **traduzir esses achados** para o paciente, sem gerar falsas certezas (estudo exploratório).

### 🧬 2.6 Enfermagem em genética/genômica
- A interpretação de painéis moleculares e biomarcadores é competência emergente da enfermagem.
- O enfermeiro atua na **comunicação de resultados**, no **aconselhamento** e na **coordenação do cuidado** baseado em perfil molecular.

---

## 3. Pontos de atenção (limitações)

| Limitação | Impacto na prática |
|-----------|--------------------|
| Estudo exploratório com dados secundários | Não usar como biomarcador validado |
| Sem validação em coorte independente | Resultados devem ser confirmados |
| Sem demonstração de causalidade | Interpretar como associação |
| Dados de expressão (RNA), não proteína | Validação funcional necessária |

> **Mensagem central:** os resultados dão **base biológica** para cuidados já
> recomendados na oncologia hepática (vigilância de TEV, controle metabólico,
> segurança medicamentosa, suporte nutricional), mas **não substituem** a
> decisão clínica nem configuram, por si, biomarcador validado.

---

## 4. Referência aos resultados do repositório

- DEGs: `results/deg/DEG_significant_LA.csv`
- Kaplan-Meier: `results/3grupos/KM_OS_*.png` e `results/3grupos/survival_logrank.csv`
- GSVA: `results/enrichment/GSVA_summary.csv`
- GSEA: `results/enrichment/GSEA_KEGG.csv`
- Reactome: `results/enrichment/Reactome_ORA_*.csv`
