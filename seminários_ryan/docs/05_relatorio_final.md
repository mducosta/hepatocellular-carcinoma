# 5. Relatório final do seminário

> Documento de síntese. Arquivos de apoio: `README.md`,
> `docs/01_analise_exploratoria.md`, `docs/02_ComBat.md`,
> `docs/03_logFC_limma.md`, `docs/04_volcano_plots.md`.

---

## 5.1 Resumo executivo

Foram analisados **529 amostras de fígado**:

- **369 tumores primários** de TCGA-LIHC;
- **50 tecidos adjacentes** de TCGA-LIHC;
- **110 fígados normais** de GTEx.

A expressão gênica foi obtida do hub Toil do UCSC Xena (log2 TPM + 0,001),
mapeada para símbolos HUGO e filtrada para **12.401 genes** expressos. O painel
da via KEGG **hsa05417** (*Lipid and Atherosclerosis*) foi analisado
separadamente (213 genes).

### Principais achados

1. **Fenótipos:** a coorte é predominantemente masculina (~67% no TCGA tumor e
   ~70% no GTEx), sem diferença significativa de sexo entre os grupos. A idade
   difere entre TCGA tumor e GTEx (59,4 vs. 53,7 anos; p = 3,96 × 10⁻⁶), o que
   motivou o ajuste por idade. A maioria dos tumores é estágio I–II (69,1%) e
   grau G2–G3 (80,8%).
2. **Estrutura global:** o UMAP mostra separação parcial entre tumor,
   adjacente e GTEx normal, e separação clara entre TCGA e GTEx. A coordenada
   UMAP2 associa-se inversamente a grau (ρ = −0,346; p = 1,2 × 10⁻¹¹) e estágio
   (ρ = −0,222; p = 3,1 × 10⁻⁵) e positivamente a idade (ρ = 0,207;
   p = 6,5 × 10⁻⁵), sugerindo que a expressão global captura um eixo de
   agressividade tumoral. UMAP não fornece variância explicada; as correlações
   são exploratórias.
3. **Expressão diferencial (FDR < 0,05 e |logFC| > 1):**
   - TCGA tumor vs. adjacente: **3.322 up / 808 down**;
   - TCGA tumor vs. adjacente (pareado): **2.934 up / 776 down**;
   - TCGA tumor vs. GTEx normal: **2.245 up / 1.143 down**;
   - TCGA adjacente vs. GTEx normal: **584 up / 1.966 down**;
   - TCGA grau G3–G4 vs. G1–G2: **61 up / 111 down**;
   - TCGA estágio III–IV vs. I–II: **3 up / 43 down**;
   - GTEx masculino vs. feminino: **8 up / 5 down**;
   - via LA (TCGA tumor vs. GTEx): **14 up / 32 down**.
4. **Biologia:** no tumor há perda de genes de identidade hepática
   (CLEC4M, CLEC4G, CLEC1B, FCN2/3, CYP1A2, GDF2, MARCO), ganho de proliferação
   (BIRC5, UBE2T, CDC25C, AURKB, CCNB1, CDC20), remodelamento vascular/matricial
   (PLVAP, MMP9, VIM, SPP1, COL1A1) e ativação de hipóxia em estágios avançados
   (CA9, MIR210HG).
5. **ComBat:** a correção de lote TCGA × GTEx é **não identificável** porque
   estudo e condição são perfeitamente confundidos. O ComBat recusa o modelo
   com covariável biológica; sem covariável, ele reduz 3.388 DEGs a 1,
   demonstrando sobrecorreção. Por isso, o ComBat **não** foi usado na análise
   principal.
6. **Controle positivo:** a comparação GTEx masculino vs. feminino recuperou
   genes do cromossomo Y e XIST, confirmando que o pipeline detecta sinal
   biológico real.

---

## 5.2 Objetivos e como foram atendidos

| Objetivo | Como foi atendido | Onde ver |
|---|---|---|
| Análise exploratória e clínica | Tabelas e figuras de sexo, idade, estágio, grau, UMAP, heatmap e marcadores | `docs/01_analise_exploratoria.md`; `resultados/tabelas/01–12`; `resultados/figuras/01–12` |
| Fenótipos avançados | Missingness, associações, AFP, fibrose, vírus, Kaplan-Meier, Cox e `cox.zph` | `docs/08_fenotipos_avancados.md`; `scripts/05_fenotipos_avancados.R`; `resultados/tabelas/13–20` |
| Auditorias, benchmarks e compliance | Auditoria de dados e do Volcano_12, benchmarks e checklist | `auditorias/` |
| TCGA vs. TCGA | Tumor vs. adjacente (não pareado e pareado); subgrupos de estágio, grau, sexo e idade | `resultados/tabelas/DE/`; `docs/04_volcano_plots.md` |
| GTEx vs. GTEx | Masculino vs. feminino; idade ≥60 vs. <60 | `resultados/tabelas/DE/`; `docs/04_volcano_plots.md` |
| TCGA vs. GTEx | Tumor vs. fígado normal; adjacente vs. GTEx normal; via LA | `resultados/tabelas/DE/`; `docs/04_volcano_plots.md` |
| Explicação do ComBat | Modelo, funcionamento, aplicação e limitação na deconvolução | `docs/02_ComBat.md`; `scripts/04_ComBat_demonstracao.R` |
| Cálculo do logFC no limma | Modelo linear, desenho, contrastes, `topTable`, exemplos | `docs/03_logFC_limma.md` |
| Volcano plots | Três principais + complementares | `resultados/figuras/Volcano_*.png`; `docs/04_volcano_plots.md` |
| Organização | Pasta `seminários_ryan/` com scripts, dados processados, resultados e docs | Este repositório |

---

## 5.3 Métodos em resumo

### 5.3.1 Dados

- **TCGA-LIHC:** expressão RNA-seq do hub Toil (UCSC Xena), fenótipos do
  GDC/TCGA e da matriz clínica `LIHC_clinicalMatrix`.
- **GTEx:** fígado normal, mesma matriz Toil; sexo e idade do GTEx v8.
- **Amostras excluídas:** 2 tumores recidivados (`-02`), por não representarem
  tumor primário.
- **Escala:** log2(TPM + 0,001).
- **Mapeamento:** Ensembl → símbolo HUGO (`org.Hs.eg.db`), isoformas agregadas
  pela média.
- **Filtro:** expressão ≥ 1 (log2 TPM) em ≥ 20% das amostras e variância > 0.

### 5.3.2 Expressão diferencial

- Pacote `limma` (Ritchie et al., 2015).
- Modelo: `~ 0 + grupo + sexo + idade` (quando aplicável).
- Comparações pareadas: `~ 0 + grupo + paciente`.
- `eBayes(trend = TRUE, robust = TRUE)`.
- Correção de múltiplos testes: Benjamini-Hochberg.
- Critério: FDR < 0,05 e |logFC| > 1.
- `logFC` em log2; fold change = 2^logFC.

### 5.3.3 ComBat

- Pacote `sva` (Johnson et al., 2007; Leek et al., 2012).
- Tentativa com `mod = ~ grupo` → erro de confundimento.
- Demonstração com `mod = NULL` → sobrecorreção.
- Não aplicado à análise principal.

### 5.3.4 Visualização

- `ggplot2`, `ggrepel`, `pheatmap`.
- Volcano plots com FDR no eixo y e logFC no eixo x.
- UMAP com os 2.000 genes mais variáveis.

---

## 5.4 Resultados detalhados

### 5.4.1 Fenótipos

| Variável | TCGA tumor | TCGA adjacente | GTEx normal |
|---|---:|---:|---:|
| n | 369 | 50 | 110 |
| Sexo feminino/masculino | 120/249 | 22/28 | 33/77 |
| Idade (média ± DP) | 59,4 ± 13,5 | 61,7 ± 16,1 | 53,7 ± 11,1 |
| Estágio I/II/III/IV | 169/86/85/5 | — | — |
| Grau G1/G2/G3/G4 | 55/177/121/11 | — | — |

- Sexo: sem diferença significativa (qui-quadrado TCGA tumor vs. GTEx:
  p = 0,703).
- Idade: diferença significativa (Wilcoxon TCGA tumor vs. GTEx:
  p = 3,96 × 10⁻⁶).
- UMAP2: grau (ρ = −0,346), estágio (ρ = −0,222), idade (ρ = 0,207).

### 5.4.2 Fenótipos avançados, sobrevida e Cox

**Missingness (TCGA tumor):** invasão vascular 70,7%, sorologia viral 55,8%,
fibrose 43,1%, Child-Pugh 35,8%, AFP 25,2%, estágio 6,5%, grau 1,4% e idade
0,3%.

**Associações com p ajustado (BH) < 0,05:**

| Associação | p | V de Cramér |
|---|---:|---:|
| Child-Pugh × estágio | 5,5 × 10⁻⁸ | 0,293 |
| log2(AFP+1) × grau (Wilcoxon) | 4,3 × 10⁻⁶ | — |
| AFP ≥400 × grau | 2,0 × 10⁻⁵ | 0,257 |
| Raça × grau | 1,0 × 10⁻⁴ | 0,251 |
| Idade × grau | 4,9 × 10⁻³ | — |
| log2(AFP+1) × sexo | 5,4 × 10⁻³ | — |
| Fibrose Ishak × estágio | 6,7 × 10⁻³ | — |
| Vírus × grau | 1,1 × 10⁻² | 0,236 |
| AFP ≥400 × estágio | 1,4 × 10⁻² | 0,152 |
| Fibrose avançada × estágio | 1,4 × 10⁻² | 0,174 |
| Vírus × estágio | 1,4 × 10⁻² | 0,233 |

Para vírus e raça, a categoria `Unknown` foi excluída dos testes.

**AFP (medianas):** estágio I–II = 12 vs. III–IV = 34,5; grau G1–G2 = 7 vs.
G3–G4 = 48; mulheres = 28 vs. homens = 11. A média é inflada por outliers;
por isso usamos mediana e log2(AFP+1).

**Sobrevida global (Kaplan-Meier):** apenas **estágio** foi significativo
(log-rank p = 5,18 × 10⁻⁷; n = 339; 114 eventos). Grau (p = 0,682), sexo
(p = 0,288), idade (p = 0,244), AFP ≥400 (p = 0,941) e fibrose avançada
(p = 0,428) não atingiram significância. Veja as figuras 17–23.

**Cox multivariado:**

| Modelo | Variável | HR (IC 95%) | p |
|---|---|---|---:|
| A (n = 337; 112 eventos) | Estágio III–IV vs. I–II | 2,61 (1,79–3,81) | 6,0 × 10⁻⁷ |
| A | Grau G3–G4 vs. G1–G2 | 1,16 (0,79–1,71) | 0,450 |
| A | Sexo masculino vs. feminino | 0,88 (0,59–1,29) | 0,508 |
| A | Idade (por ano) | 1,01 (0,996–1,03) | 0,155 |
| B (n = 179; 52 eventos) | Estágio III–IV vs. I–II | 2,23 (1,24–4,03) | 7,7 × 10⁻³ |
| B | Idade (por ano) | 1,03 (1,006–1,058) | 0,0165 |
| B | AFP ≥400 vs. <400 | 1,04 (0,51–2,12) | 0,907 |
| B | Fibrose avançada vs. não | 1,16 (0,62–2,17) | 0,645 |

A suposição de riscos proporcionais foi globalmente satisfeita
(`cox.zph`: p = 0,386 no Modelo A e p = 0,153 no Modelo B); AFP no Modelo B
ficou no limite (p = 0,049).

**Correlações de Spearman:** AFP × grau (ρ = 0,292; p = 8,6 × 10⁻⁷),
AFP × estágio (ρ = 0,162; p = 0,0083), fibrose × estágio (ρ = −0,205;
p = 0,0038), idade × grau (ρ = −0,131; p = 0,0126) e idade × estágio
(ρ = −0,122; p = 0,0235).

**Interpretação:** o estágio é o fator prognóstico dominante; o grau perde
significância após ajuste pelo estágio. AFP se associa mais a grau do que a
estágio. A missingness elevada de AFP, fibrose, sorologia e invasão vascular
reduz o poder dos modelos B e das curvas correspondentes; esses resultados são
exploratórios.

### 5.4.3 DEGs

Tabela completa em `resultados/tabelas/DE/RESUMO_DEGs_TODAS_COMPARACOES.csv`.

| Comparação | Up | Down | Total |
|---|---:|---:|---:|
| TCGA tumor vs. adjacente | 3.322 | 808 | 4.130 |
| TCGA tumor vs. adjacente (pareado) | 2.934 | 776 | 3.710 |
| TCGA tumor vs. GTEx normal | 2.245 | 1.143 | 3.388 |
| TCGA adjacente vs. GTEx normal | 584 | 1.966 | 2.550 |
| TCGA masculino vs. feminino | 155 | 64 | 219 |
| TCGA G3–G4 vs. G1–G2 | 61 | 111 | 172 |
| GTEx ≥60 vs. <60 | 2 | 113 | 115 |
| TCGA ≥65 vs. <65 | 43 | 30 | 73 |
| TCGA estágio III–IV vs. I–II | 3 | 43 | 46 |
| Via LA: TCGA tumor vs. GTEx | 14 | 32 | 46 |
| Via LA: TCGA tumor vs. adjacente | 11 | 33 | 44 |
| GTEx masculino vs. feminino | 8 | 5 | 13 |

### 5.4.4 Interpretação biológica integrada

1. **Desdiferenciação hepática:** queda de CLEC4M, CLEC4G, CLEC1B, FCN2,
   FCN3, CYP1A2, GDF2, MARCO, CRHBP e ADAMTS13 no tumor. Isso é consistente com
   perda de fenótipo hepatocitário e alteração das células endoteliais
   sinusoidais e de Kupffer.
2. **Proliferação:** aumento de BIRC5, UBE2T, CDC25C, AURKB, CCNB1, CDC20,
   NUF2, KIF2C e CDKN3. O grau histológico captura bem esse eixo.
3. **Remodelamento vascular/matricial:** aumento de PLVAP, CD34, MMP9, MMP1,
   VIM, SPP1, COL1A1 e VCAM1.
4. **Hipóxia e estágio avançado:** CA9, MIR210HG e SLC6A8 sobem em III–IV;
   genes metabólicos (GNMT, CYP2C9, ACOT12, FTCD, HPX, SLC10A1) caem.
5. **Via LA:** ativação de componentes inflamatórios/oxidativos (PYCARD, LY96,
   NOX1, VCAM1, MMP1/MMP9) e alteração de receptores lipídicos (CD36, ABCG1,
   PPARG).
6. **Controle de sexo:** genes Y e XIST confirmam a qualidade dos metadados.

---

### 5.4.5 Auditorias, benchmarks e compliance

- **Auditoria de dados:** 22 verificações; todas OK, exceto o alerta de
  confundimento estrutural estudo × condição (Cramér's V = 1). Nenhum NA,
  infinito, gene duplicado ou amostra dessincronizada na matriz. Ver
  `auditorias/AUDITORIA_DADOS.md`.
- **Auditoria do Volcano_12:** a figura original não estava tecnicamente em
  branco (≈3,65% de pixels não brancos), mas apenas 1 DEG permaneceu após o
  ComBat. A versão corrigida usa dois painéis (antes/depois), tem ≈9,48% de
  pixels não brancos e está documentada em
  `auditorias/AUDITORIA_VOLCANO_12.md`.
- **Benchmarks:** carga da matriz em ~1,5 s; UMAP em ~4,3 s; limma com
  covariáveis em ~0,7 s (12.401 genes × 478 amostras). Concordância:
  com vs. sem covariáveis Jaccard = 0,948 e ρ = 0,999; pareado vs. não
  pareado Jaccard = 0,764 e ρ = 0,969; ComBat antes vs. depois Jaccard = 0 e
  ρ = 0,647. Ver `auditorias/BENCHMARKS.md`.
- **Compliance:** checagem de arquivos, PNGs não vazios, ausência de PDFs e
  de métodos lineares de redução de dimensionalidade, CSVs, scripts parseáveis, sementes, links e metadados. Ver
  `auditorias/COMPLIANCE.md`.

## 5.5 Limitações e ameaças à validade

1. **Confundimento TCGA × GTEx.** Estudo e condição são perfeitamente
   confundidos; o `logFC` tumor vs. normal combina biologia e lote. A
   comparação TCGA adjacente vs. GTEx normal mostra que os tecidos normais já
   diferem entre estudos (2.550 DEGs), reforçando a cautela.
2. **Composição celular.** O bulk mede a média de tipos celulares; parte das
   diferenças reflete proporções celulares, não expressão intrínseca. A
   deconvolução celular não foi realizada neste seminário; é uma extensão
   natural.
3. **TPM e composicionalidade.** TPM é relativo e não preserva contagens
   absolutas; o `limma` em log2 TPM é uma aproximação. Para contagens, `voom`
   ou modelos binomial negativo seriam mais adequados.
4. **Filtro de genes.** O filtro de expressão e o mapeamento Ensembl → símbolo
   removem genes; os resultados dependem dessas escolhas.
5. **Poder estatístico.** Subgrupos como estágio IV (n = 5) e grau G4 (n = 11)
   têm poder limitado; os resultados de subgrupos são exploratórios.
6. **Idade no GTEx.** A idade é uma faixa convertida para o ponto médio, o que
   introduz erro de medida.
7. **Ausência de validação experimental.** Não houve qPCR, IHC ou coorte
   independente; os achados são computacionais.
8. **ComBat não é solução para confundimento.** A demonstração mostra que a
   correção cega pode eliminar a biologia.
9. **Missingness clínico.** Invasão vascular (70,7%), sorologia viral (55,8%),
   fibrose (43,1%), Child-Pugh (35,8%) e AFP (25,2%) têm dados ausentes ou
   desconhecidos; os modelos que usam essas variáveis têm poder limitado.

---

## 5.6 Conclusões

- O pipeline integra TCGA-LIHC e GTEx de forma reprodutível e documentada,
  com controle de qualidade, análise clínica, expressão diferencial e
  visualização.
- A comparação **TCGA tumor vs. adjacente** é a evidência interna mais sólida
  de reprogramação transcricional no HCC: perda de identidade hepática, ganho
  de proliferação e remodelamento vascular/matricial.
- A comparação **TCGA vs. GTEx** confirma e amplifica esses achados, mas deve
  ser interpretada com cautela devido ao confundimento com o estudo.
- A comparação **GTEx masculino vs. feminino** funciona como controle positivo
  e valida a anotação de sexo e o pipeline.
- O **ComBat** é uma ferramenta útil quando lote e biologia são separáveis;
  neste conjunto, a confusão estrutural impede sua aplicação à comparação
  TCGA × GTEx. Na deconvolução celular, seu uso exige incluir a composição
  celular no modelo e validar com marcadores.
- O **logFC do limma** é uma diferença de médias em escala log2, obtida por
  modelo linear e contraste; sua interpretação correta depende do contraste,
  das covariáveis e da escala dos dados.
- Nos fenótipos avançados, **estágio** foi o fator prognóstico dominante
  (log-rank p = 5,18 × 10⁻⁷; Cox multivariado HR = 2,61), enquanto **grau**
  perdeu significância após ajuste. **AFP** associou-se a grau, e **Child-Pugh**
  e **fibrose** associaram-se a estágio, refletindo doença hepática de base.
- As auditorias confirmaram a integridade dos dados e dos artefatos; o único
  alerta foi o confundimento estrutural TCGA × GTEx. O `Volcano_12` foi
  auditado e refeito como painel antes/depois.

---

## 5.7 Recomendações para próximos passos

1. **Deconvolução celular** com referência de fígado (ex.: CIBERSORTx, MuSiC,
   BayesPrism), seguida de associação das proporções com estágio/grau.
2. **Validação em coorte independente** (ex.: ICGC LIRI-JP, GSE14520) para os
   principais DEGs.
3. **Análise de sobrevida** dos genes candidatos (OS, PFI) com modelos de Cox
   ajustados por estágio, grau, sexo, idade e função hepática; incorporar
   Child-Pugh, fibrose e AFP em coortes com menos missingness.
4. **Enriquecimento funcional** (GSEA/GSVA) para interpretar os módulos de
   expressão associados a grau e estágio.
5. **Correção de lote com desenho adequado** (ex.: ComBat dentro de cada
   estudo, ou métodos que preservem composição celular) em vez de corrigir
   TCGA × GTEx diretamente.
6. **Análise de células imunes** (ssGSEA, xCell, CIBERSORTx) para investigar
   o microambiente do HCC.

---

## 5.8 Referências principais

- Ritchie, M. E. et al. **limma powers differential expression analyses for
  RNA-sequencing and microarray studies.** *Nucleic Acids Research*, 2015.
- Johnson, W. E.; Li, C.; Rabinovic, A. **Adjusting batch effects in microarray
  expression data using empirical Bayes methods.** *Biostatistics*, 2007.
- Leek, J. T. et al. **The sva package for removing batch effects and other
  unwanted variation in high-throughput experiments.** *Bioinformatics*, 2012.
- Zhang, Y.; Parmigiani, G.; Johnson, W. E. **ComBat-seq: batch effect
  adjustment for RNA-seq count data.** *NAR Genomics and Bioinformatics*, 2020.
- Goldman, M. J. et al. **Visualizing and interpreting cancer genomics data via
  the Xena platform.** *Nature Biotechnology*, 2020.
- GTEx Consortium. **The GTEx Consortium atlas of genetic regulatory effects
  across human tissues.** *Science*, 2020.
- Cancer Genome Atlas Research Network. **Comprehensive and integrative
  genomic characterization of hepatocellular carcinoma.** *Cell*, 2017.
- Newman, A. M. et al. **Determining cell type abundance and expression from
  bulk tissues with digital cytometry.** *Nature Biotechnology*, 2019.
- Wang, X. et al. **Bulk tissue cell type deconvolution with multi-subject
  single-cell expression reference.** *Nature Communications*, 2019.
- Korsunsky, I. et al. **Fast, sensitive and accurate integration of
  single-cell data with Harmony.** *Nature Methods*, 2019.
