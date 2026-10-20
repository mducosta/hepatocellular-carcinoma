# Relatório consolidado — Seminário Ryan (TCGA-LIHC × GTEx)

> Documento gerado a partir de `docs/*.md`.
> Versão **somente Markdown**; não há HTML consolidado.
> Para regenerar: `python scripts/09_gerar_relatorio_md.py`.

## Sumário

- [01_analise_exploratoria](docs/01_analise_exploratoria.md)
- [02_ComBat](docs/02_ComBat.md)
- [03_logFC_limma](docs/03_logFC_limma.md)
- [04_volcano_plots](docs/04_volcano_plots.md)
- [05_relatorio_final](docs/05_relatorio_final.md)
- [06_dicionario_dados](docs/06_dicionario_dados.md)
- [07_reproducibilidade](docs/07_reproducibilidade.md)
- [08_fenotipos_avancados](docs/08_fenotipos_avancados.md)
- [09_auditorias_benchmarks](docs/09_auditorias_benchmarks.md)
- [10_grupos_alternativos](docs/10_grupos_alternativos.md)
- [11_grupos_separados](docs/11_grupos_separados.md)


---

# 1. Análise exploratória e fenótipos clínicos

> Arquivo: `docs/01_analise_exploratoria.md`.
> Script: `scripts/01_analise_exploratoria.R`.
> Saídas: `resultados/tabelas/01–12_*.csv` e `resultados/figuras/01–12_*.png`.

---

## 1.1 Dados utilizados

| Item | TCGA-LIHC | GTEx (fígado) |
|---|---|---|
| Tumores primários | 369 | — |
| Tecido adjacente (solid tissue normal) | 50 | — |
| Fígado normal | — | 110 |
| Tumores recidivados | 2 (excluídos das análises principais) | — |
| **Total analisado** | **419** | **110** |

- **Expressão:** hub Toil do UCSC Xena, `TcgaTargetGtex_rsem_gene_tpm`,
  em log2(TPM + 0,001).
- **Mapeamento:** IDs Ensembl → símbolo HUGO (`org.Hs.eg.db`); isoformas
  agregadas pela média.
- **Filtro:** gene detectado (log2 TPM ≥ 1) em ≥ 20% das amostras e variância
  > 0.
- **Matriz final:** 12.401 genes × 529 amostras (369 tumor, 50 adjacente,
  110 GTEx normal).
- **Painel KEGG hsa05417** (*Lipid and Atherosclerosis*): 213 genes,
  analisado separadamente.

### Fontes clínicas

- TCGA-LIHC: `LIHC_clinicalMatrix` (UCSC Xena/TCGA) e
  `TCGA-LIHC.clinical.tsv.gz` (GDC).
- GTEx: `GTEX_phenotype` e
  `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` (GTEx v8).
- Variáveis: sexo, idade, estágio AJCC, grau histológico, tipo de amostra,
  status vital, AFP, fibrose de Ishak e sorologia viral.

Os metadados consolidados estão em:

- `dados/processed/metadata_amostras.tsv` (todas as 531 amostras, incluindo
  recidivas);
- `dados/processed/metadata_analise.tsv` (529 amostras, sem recidivas).

---

## 1.2 Distribuição das amostras

| Grupo | n |
|---|---:|
| TCGA — tumor primário | 369 |
| TCGA — tecido adjacente | 50 |
| GTEx — fígado normal | 110 |

![Amostras por grupo](resultados/figuras/01_amostras_por_grupo.png)

*Figura 1. Distribuição das amostras hepáticas por grupo (apenas grupos com n > 0).*

> A versão anterior cruzava estudo × condição e exibia combinações vazias
> (`GTEx × Tumor = 0`, `TCGA × Normal = 0` etc.). A figura atual conta apenas
> os grupos realmente observados; as 2 recidivas, excluídas da análise
> principal, aparecem na Figura 26 (`docs/10_grupos_alternativos.md`).

---

## 1.3 Fenótipos clínicos

### 1.3.1 Sexo

| Grupo | Feminino | Masculino | Total |
|---|---:|---:|---:|
| TCGA tumor | 120 (32,5%) | 249 (67,5%) | 369 |
| TCGA adjacente | 22 (44,0%) | 28 (56,0%) | 50 |
| GTEx normal | 33 (30,0%) | 77 (70,0%) | 110 |

Testes qui-quadrado:

| Comparação | p |
|---|---:|
| TCGA tumor vs. GTEx normal | 0,703 |
| TCGA tumor vs. TCGA adjacente | 0,147 |
| TCGA adjacente vs. GTEx normal | 0,121 |

**Interpretação:** a proporção de homens é alta em todos os grupos, mas não há
diferença estatisticamente significativa entre eles. Ainda assim, sexo foi
incluído como covariável nas análises de expressão diferencial, pois pode
explicar parte da variância (especialmente genes ligados ao cromossomo X/Y).

![Sexo por grupo](resultados/figuras/02_sexo_por_grupo.png)

*Figura 2. Composição por sexo em cada grupo.*

### 1.3.2 Idade

| Grupo | n | Média ± DP | Mediana | Q1–Q3 | Mín–Máx |
|---|---:|---:|---:|---:|---:|
| TCGA tumor | 368 | 59,4 ± 13,5 | 61 | 51,8–69,0 | 16–90 |
| TCGA adjacente | 50 | 61,7 ± 16,1 | 68,5 | 52,3–74,5 | 20–81 |
| GTEx normal | 110 | 53,7 ± 11,1 | 55 | 45,0–65,0 | 25–65* |

\* A idade no GTEx é derivada do ponto médio da faixa etária do doador
(`20-29`, `30-39`, …, `60-69`), por isso os valores são discretizados.

Testes:

| Comparação | Teste | p |
|---|---|---:|
| TCGA tumor vs. GTEx normal | Wilcoxon | 3,96 × 10⁻⁶ |
| TCGA tumor vs. TCGA adjacente | Wilcoxon | 0,0438 |
| Três grupos | Kruskal-Wallis | 8,80 × 10⁻⁷ |

**Interpretação:** os tumores TCGA são, em média, ~5,7 anos mais velhos que os
doadores GTEx. Essa diferença é estatisticamente significativa e justifica o
ajuste por idade nas comparações TCGA × GTEx. A idade também foi usada para
comparações internas (≥65 vs. <65 no TCGA; ≥60 vs. <60 no GTEx).

![Idade por grupo](resultados/figuras/03_idade_por_grupo.png)

*Figura 3. Distribuição de idade por grupo (violino + boxplot).*

### 1.3.3 Estágio AJCC (TCGA tumor)

| Estágio | n | % |
|---|---:|---:|
| I | 169 | 45,8 |
| II | 86 | 23,3 |
| III | 85 | 23,0 |
| IV | 5 | 1,4 |
| Não informado | 24 | 6,5 |

![Estágio TCGA](resultados/figuras/04_estagio_TCGA.png)

*Figura 4. Estágio AJCC dos tumores primários TCGA-LIHC.*

**Interpretação:** há predomínio de estágios iniciais (I–II = 69,1%), o que é
típico de coortes cirúrgicas. O número de tumores estágio IV é pequeno (n = 5),
limitando a comparação IV vs. I–III; por isso a análise de expressão agrupou
III–IV vs. I–II.

### 1.3.4 Grau histológico (TCGA tumor)

| Grau | n | % |
|---|---:|---:|
| G1 | 55 | 14,9 |
| G2 | 177 | 48,0 |
| G3 | 121 | 32,8 |
| G4 | 11 | 3,0 |
| Não informado | 5 | 1,4 |

![Grau TCGA](resultados/figuras/05_grau_TCGA.png)

*Figura 5. Grau histológico dos tumores primários TCGA-LIHC.*

**Interpretação:** a maioria é G2–G3. Para aumentar o poder, a análise de
expressão comparou G3–G4 vs. G1–G2.

### 1.3.5 Estágio × grau

![Estágio x grau](resultados/figuras/06_estagio_x_grau_TCGA.png)

*Figura 6. Distribuição conjunta de estágio agrupado e grau agrupado.*

> A exploração ampliada de fenótipos — missingness, AFP, fibrose, sorologia
> viral, Child-Pugh, histologia, invasão vascular, Kaplan-Meier e Cox — está em
> [`docs/08_fenotipos_avancados.md`](docs/08_fenotipos_avancados.md).

---

## 1.4 Padrões globais de expressão (UMAP)

Foi feito **UMAP** (*Uniform Manifold Approximation and Projection*) com os
**2.000 genes mais variáveis**, usando:

| Parâmetro | Valor |
|---|---|
| Genes | 2.000 mais variáveis |
| `n_neighbors` | 30 |
| `min_dist` | 0,25 |
| Métrica | cosine |
| `scale` | TRUE |
| `pca` (parâmetro interno do `uwot`) | 50 |
| Semente | 2026 |

Diferentemente de métodos lineares de redução de dimensionalidade, o UMAP **não fornece percentual de variância
explicada** e seus eixos não têm interpretação linear direta. Ele é usado
como visualização de vizinhança e estrutura local, não como coordenada
quantitativa.

![UMAP por grupo](resultados/figuras/07_UMAP_grupo.png)

*Figura 7. UMAP colorido por grupo. TCGA tumor, TCGA adjacente e GTEx normal
separam-se parcialmente.*

![UMAP por estudo](resultados/figuras/08_UMAP_estudo.png)

*Figura 8. UMAP colorido por estudo. A separação TCGA vs. GTEx é evidente,
indicando forte efeito de estudo/lote.*

![UMAP por sexo (TCGA tumor)](resultados/figuras/09_UMAP_sexo_TCGA.png)

*Figura 9. UMAP dentro dos tumores TCGA, colorido por sexo.*

![UMAP estágio e grau (TCGA tumor)](resultados/figuras/10_UMAP_estagio_grau_TCGA.png)

*Figura 10. UMAP dos tumores TCGA, colorido por estágio e com forma por grau.*

### Associação das coordenadas UMAP com a clínica

Correlações de Spearman dentro dos tumores TCGA (apenas exploratórias):

| Variável | UMAP1 (ρ) | p (UMAP1) | UMAP2 (ρ) | p (UMAP2) |
|---|---:|---:|---:|---:|
| Idade | −0,140 | 0,0073 | 0,207 | 6,5 × 10⁻⁵ |
| Estágio | 0,076 | 0,159 | −0,222 | 3,1 × 10⁻⁵ |
| Grau | 0,108 | 0,040 | −0,346 | 1,2 × 10⁻¹¹ |

**Interpretação:**

- **UMAP2** captura um eixo inversamente associado à **agressividade tumoral**:
  quanto maior o grau, menor a coordenada UMAP2 (ρ = −0,346); o mesmo padrão,
  mais fraco, aparece para estágio (ρ = −0,222).
- **UMAP2** também se associa positivamente com idade (ρ = 0,207), indicando
  que idade e biologia tumoral contribuem de formas parcialmente
  independentes.
- **UMAP1** tem associações fracas (idade ρ = −0,140; grau ρ = 0,108; estágio
  ρ = 0,076) e não deve ser interpretado como um eixo clínico.
- Como o UMAP é não linear, essas correlações servem apenas como pistas
  exploratórias; a inferência principal vem do `limma` e dos modelos de
  sobrevida.

---

## 1.5 Heatmap dos genes mais variáveis

![Heatmap top 50](resultados/figuras/11_heatmap_top50_variaveis.png)

*Figura 11. Heatmap dos 50 genes mais variáveis (z-score por gene). As colunas
são amostras; as barras superiores indicam grupo e sexo.*

**Interpretação:** há blocos de genes com expressão claramente distinta entre
TCGA tumor e GTEx normal, e também heterogeneidade dentro dos tumores. A
estrutura por grupo é consistente com o UMAP e com os resultados de expressão
diferencial.

---

## 1.6 Marcadores conhecidos de HCC e remodelamento tecidual

![Violinos de marcadores](resultados/figuras/12_violino_marcadores.png)

*Figura 12. Expressão de marcadores selecionados nos três grupos (violino + boxplot).*

Destaques:

- **GPC3** e **AFP**: marcadores clássicos de HCC, com expressão elevada no
  tumor.
- **BIRC5** (survivina): gene de proliferação/antiapoptose, elevado no tumor.
- **MMP9**, **COL1A1**, **VIM**, **SPP1**: associados a remodelamento da matriz
  extracelular, invasão e transição epitélio-mesenquimal.
- **CDH1** (E-caderina): tende a ser menor no tumor, consistente com perda de
  adesão epitelial.
- **KRT19** e **EPCAM**: marcadores de progenitores/ductais; a interpretação é
  mais heterogênea e depende do subtítulo molecular.

Os valores individuais estão em
`resultados/tabelas/12_expressao_marcadores.csv`.

---

## 1.7 Síntese da análise exploratória

1. **A coorte é predominantemente masculina** (≈67% no TCGA tumor e ≈70% no
   GTEx), sem diferença significativa entre os grupos.
2. **A idade difere entre TCGA e GTEx** (p = 3,96 × 10⁻⁶), o que motivou o
   ajuste por idade nas comparações entre estudos.
3. **A maioria dos tumores é estágio I–II e grau G2–G3**; os subgrupos extremos
   (estágio IV, grau G4) são pequenos.
4. **UMAP2 resume um eixo de agressividade** (grau > estágio > idade) e mostra
   que a expressão global carrega informação clínica relevante.
5. **Há separação global entre TCGA e GTEx** (UMAP), compatível com efeito de
   estudo/lote e com diferenças de composição tecidual.
6. **Marcadores clássicos de HCC** (GPC3, AFP, BIRC5) e de remodelamento
   (MMP9, COL1A1, VIM, SPP1) apresentam o comportamento esperado, validando o
   processamento e a anotação.

---

## 1.8 Limitações desta etapa

- A idade do GTEx é uma **faixa** convertida para o ponto médio; isso introduz
  erro de medida e reduz a precisão das comparações etárias.
- O estágio e o grau são **do paciente**, não necessariamente da amostra
  específica; no TCGA, os 50 tecidos adjacentes herdam o estágio/grau do tumor
  do mesmo paciente. Por isso as análises de estágio/grau usam apenas tumores.
- O UMAP é sensível à seleção de genes (top 2.000 por variância), aos
  hiperparâmetros (`n_neighbors`, `min_dist`) e ao efeito de lote. A separação
  por estudo não deve ser interpretada como prova de diferença biológica.
- O UMAP não fornece variância explicada e seus eixos não são lineares; por
  isso, as correlações com clínica são apenas exploratórias.
- A análise exploratória é **descritiva**; os testes de hipótese principais
  estão na etapa de expressão diferencial (`docs/04_volcano_plots.md`).

---

# 2. ComBat (`sva`): o que é, como funciona e limites na deconvolução celular

> Documento de apoio ao seminário. Arquivo: `docs/02_ComBat.md`.
> Demonstração executável: `scripts/04_ComBat_demonstracao.R`.
> Resultados: `resultados/tabelas/DE/COMPARACAO_ComBat_antes_depois.csv`,
> `resultados/qc/resumo_ComBat_demonstracao.txt`.

---

## 2.1 O que é o ComBat?

O **ComBat** é uma função do pacote Bioconductor **`sva`** (*Surrogate Variable
Analysis*) voltada à **correção de efeitos de lote** (*batch effects*) em dados
de expressão gênica de alta dimensão. O nome vem de *Combating Batch Effects*
(“combatendo efeitos de lote”).

Um **efeito de lote** é uma fonte de variação técnica ou de amostragem que não
interessa biologicamente, mas que altera a expressão medida. Exemplos:

- amostras processadas em dias, placas, centros ou sequenciadores diferentes;
- coortes distintas (por exemplo, TCGA e GTEx);
- diferentes protocolos de extração de RNA, bibliotecas ou pipelines;
- diferentes technologias (microarray vs. RNA-seq, RNA-seq bulk vs. single-cell).

Se não for tratado, o lote pode:

1. criar agrupamentos artificiais em análises de redução de dimensionalidade (ex.: UMAP);
2. inflar ou reduzir taxas de falsos positivos na expressão diferencial;
3. fazer com que a variável de interesse fique confundida com a variável
   técnica, tornando a interpretação inválida.

O ComBat tenta remover o efeito de lote **preservando** as diferenças
biológicas que forem explicitamente informadas no modelo.

---

## 2.2 Modelo estatístico

Para cada gene $g$ e cada amostra $i$, o ComBat assume um modelo aditivo e
multiplicativo:

$$
y_{gi} = \alpha_g + X_i\beta_g + \gamma_{g,b(i)} + \delta_{g,b(i)}\,\varepsilon_{gi}
$$

em que:

| Termo | Interpretação |
|---|---|
| $y_{gi}$ | expressão do gene $g$ na amostra $i$ (em escala logarítmica) |
| $\alpha_g$ | expressão média global do gene $g$ |
| $X_i$ | vetor de covariáveis biológicas da amostra $i$ (ex.: grupo tumoral, sexo, idade) |
| $\beta_g$ | efeitos das covariáveis biológicas no gene $g$ |
| $b(i)$ | lote ao qual a amostra $i$ pertence |
| $\gamma_{g,b(i)}$ | efeito **aditivo** (deslocamento da média) do lote $b$ no gene $g$ |
| $\delta_{g,b(i)}$ | efeito **multiplicativo** (escala/dispersão) do lote $b$ no gene $g$ |
| $\varepsilon_{gi}$ | erro residual, em geral assumido normal com média zero |

O ComBat estima, para cada gene e cada lote, os parâmetros de **locação**
($\gamma$) e **escala** ($\delta$). Em seguida, os valores são ajustados para
voltar à escala comum:

$$
y^{*}_{gi} = \frac{y_{gi} - \hat{\alpha}_g - X_i\hat{\beta}_g - \hat{\gamma}_{g,b(i)}}{\hat{\delta}_{g,b(i)}} + \hat{\alpha}_g + X_i\hat{\beta}_g
$$

Na prática, o ComBat implementa uma versão **empírico-bayesiana**: os
parâmetros de cada gene não são estimados de forma totalmente independente.
Assume-se uma distribuição *a priori* comum para os efeitos de lote e a
estimativa de cada gene é “encolhida” (*shrinkage*) em direção a essa
distribuição. Isso estabiliza as estimativas, sobretudo quando o número de
amostras por lote é pequeno.

### 2.2.1 Principais argumentos da função

```r
sva::ComBat(
  dat       = matriz_log2,       # genes x amostras
  batch     = vetor_de_lote,     # fator com o lote de cada amostra
  mod       = model.matrix(~ grupo + sexo + idade, data = metadados),
  par.prior = TRUE,              # estimação paramétrica empírico-bayesiana
  mean.only = FALSE,             # TRUE = corrige apenas a média (locação)
  ref.batch = NULL               # opcional: escolhe um lote de referência
)
```

- **`dat`**: matriz numérica em escala logarítmica (log2 de TPM, logCPM,
  intensidades normalizadas etc.). O ComBat **não** foi desenhado para
  contagens brutas.
- **`batch`**: fator ou vetor que identifica o lote.
- **`mod`**: matriz de desenho com as variáveis biológicas que devem ser
  **preservadas**. É o argumento mais importante para a validade da análise.
- **`par.prior = TRUE`**: usa prioris paramétricas (mais comum e mais rápido).
  Com `FALSE`, usa estimação não paramétrica.
- **`mean.only = TRUE`**: corrige somente o deslocamento médio, sem mexer na
  variância. Útil quando se quer uma correção mais conservadora.
- **`ref.batch`**: define um lote de referência; os demais são ajustados a ele.
  Útil quando um lote é claramente o “padrão-ouro”.

Existe também o **`ComBat-seq`** (pacote `sva`), apropriado para **contagens
brutas** de RNA-seq, que modela os dados com distribuição binomial negativa.
Para TPM/log2 TPM, o ComBat clássico é o mais usado.

---

## 2.3 O que o ComBat **não** faz

É fundamental entender as limitações:

1. **Não resolve confundimento estrutural.** Se todos os tumores vêm do TCGA e
   todos os normais vêm do GTEx, então “estudo” e “condição tumor/normal” são
   perfeitamente confundidos. O ComBat não consegue separar o efeito de lote do
   efeito biológico, porque não há informação suficiente para isso.
2. **Não substitui um bom desenho experimental.** Nenhum método estatístico
   recupera a biologia perdida por um desenho totalmente confundido.
3. **Pode remover sinal biológico** se a composição celular ou a condição
   biológica não estiverem representadas em `mod`.
4. **Não é indicado para contagens brutas.** Use `ComBat-seq` ou métodos
   específicos para contagens.
5. **Pode distorcer genes com efeito de lote não linear** ou genes com
   expressão muito baixa/zero.
6. **É sensível ao número de amostras por lote.** Lotes muito pequenos geram
   estimativas instáveis, mesmo com *shrinkage*.

---

## 2.4 Demonstração com os dados TCGA-LIHC × GTEx

A demonstração está em `scripts/04_ComBat_demonstracao.R`. Foram usadas 478
amostras (368 tumores TCGA e 110 fígados normais GTEx) e 12.401 genes.

### Cenário 1 — tentativa de preservar a condição

Foi solicitado:

```r
mod <- model.matrix(~ group_tg, data = meta_cb)  # group_tg = tumor/normal
ComBat(dat = X, batch = study, mod = mod)
```

Resultado:

```text
Error: The covariate is confounded with batch! Remove the covariate and rerun ComBat
```

Ou seja, o próprio ComBat detecta que a condição biológica (tumor vs. normal)
é **perfeitamente confundida** com o lote (TCGA vs. GTEx). Nesse cenário, não
existe informação suficiente para estimar simultaneamente o efeito de lote e o
efeito biológico.

### Cenário 2 — ComBat sem covariável biológica

Para fins didáticos, rodamos o ComBat **sem** informar a condição:

```r
X_cor <- ComBat(dat = X, batch = study, mod = NULL)
```

O lote é removido, mas como a condição é indistinguível do lote, a diferença
tumor-normal também é removida. Resultados:

| Métrica | Antes do ComBat | Depois do ComBat (sem covariável) |
|---|---:|---:|
| DEGs (FDR < 0,05 e \|logFC\| > 1) | 3.388 | **1** |
| Up-regulados | 2.245 | 1 |
| Down-regulados | 1.143 | 0 |
| DEGs em comum | — | 1 |
| Jaccard dos DEGs | — | 0,000 |
| Correlação de logFC (Spearman) | — | 0,647 |

**Interpretação:** o ComBat “funcionou” tecnicamente — removeu o lote —, mas
removeu também a biologia. O resultado (1 DEG) é um exemplo de
**sobrecorreção** (*overcorrection*). Isso ilustra por que não aplicamos ComBat
cegamente à comparação TCGA × GTEx.

![UMAP antes e depois do ComBat](resultados/figuras/13_UMAP_ComBat_antes_depois.png)

*Figura 13. UMAP com os 2.000 genes mais variáveis antes e depois do ComBat sem
covariável. A separação entre TCGA e GTEx desaparece, mas a diferença
tumor-normal também.*

![AFP antes e depois do ComBat](resultados/figuras/14_AFP_ComBat_antes_depois.png)

*Figura 14. Exemplo com AFP (violino + boxplot): a expressão média cai e a
diferença biológica esperada entre tumor e fígado normal é atenuada após a
correção sem covariável.*

### O que fazemos na análise principal

Na análise principal (`scripts/02_expressao_diferencial.R`), **não aplicamos
ComBat** à comparação TCGA tumor × GTEx normal. Em vez disso:

- usamos `limma` com **sexo e idade como covariáveis**;
- tratamos a comparação **TCGA tumor × TCGA adjacente** como a comparação
  interna mais limpa para identificar genes associados à transformação
  maligna;
- usamos **TCGA adjacente × GTEx normal** para estimar o quanto os tecidos
  “normais” diferem entre os dois estudos (uma medida do efeito de lote e de
  diferenças de composição tecidual);
- interpretamos a comparação TCGA × GTEx com cautela, sabendo que o efeito
  estimado combina biologia tumoral + diferenças de estudo/tecido.

---

## 2.5 Aplicação e limitação do ComBat na deconvolução celular

### 2.5.1 O problema da deconvolução

Dados de expressão **bulk** (como TCGA e GTEx) medem a média de milhares de
células. A **deconvolução celular** tenta estimar a proporção de cada tipo
celular (hepatócitos, células imunes, fibroblastos, células endoteliais etc.)
a partir dessa mistura, usando assinaturas de referência (scRNA-seq ou sorted
cells). Exemplos de métodos: CIBERSORTx, MuSiC, BayesPrism, cell2location,
Tangram, entre outros.

O problema é que a deconvolução é extremamente sensível a:

- diferenças de plataforma entre a referência e o alvo (microarray vs. RNA-seq,
  bulk vs. single-cell);
- diferenças de protocolo e de estudo;
- genes marcadores de tipo celular que também são afetados por lote;
- composição celular diferente entre lotes.

### 2.5.2 Onde o ComBat pode ajudar

O ComBat pode ser aplicado **antes** da deconvolução para harmonizar:

- matrizes de referência construídas em estudos diferentes;
- amostras bulk de coortes diferentes;
- pseudobulk de scRNA-seq de diferentes doadores/batches.

Se a composição celular for conhecida ou puder ser aproximada, ela pode ser
incluída em `mod` para que o ComBat **preserve** as diferenças de composição.
Também é possível corrigir apenas a média (`mean.only = TRUE`) para reduzir o
risco de distorcer a variabilidade biológica.

### 2.5.3 Por que o ComBat pode ser perigoso na deconvolução

O ComBat assume que o efeito de lote é **independente** da biologia. Na
deconvolução, essa premissa frequentemente falha:

1. **A composição celular varia com o lote.** Por exemplo, amostras processadas
   em um centro podem ter mais infiltrado imune; se isso não estiver em `mod`,
   o ComBat pode remover justamente o sinal de composição que se quer estimar.
2. **Assinaturas de tipos celulares podem ser confundidas com lote.** Genes
   marcadores de hepatócitos, células T, macrófagos etc. podem ter efeitos de
   lote específicos. Corrigi-los pode apagar a assinatura que a deconvolução
   precisa.
3. **Lote e condição podem ser confundidos.** Como no exemplo TCGA × GTEx, se
   todos os tumores vêm de um estudo e todos os normais de outro, o ComBat não
   consegue separar lote de doença.
4. **Correção em escala logarítmica não é neutra.** A transformação
   aditiva/multiplicativa pode alterar a relação entre genes e a composição
   celular, sobretudo para genes de baixa expressão.
5. **Referência e alvo podem ter distribuições diferentes.** O ComBat pode
   aproximar as distribuições, mas não garante que a relação linear assumida
   pela deconvolução seja preservada.

### 2.5.4 Recomendações práticas

- **Inclua a composição celular em `mod`** quando ela for conhecida ou puder ser
  estimada (por exemplo, com um método robusto inicial).
- **Faça a correção por lote dentro de cada condição biológica**, quando
  possível, em vez de misturar condições confundidas.
- **Use `mean.only = TRUE`** se o objetivo for apenas alinhar médias e preservar
  a variabilidade.
- **Valide com marcadores conhecidos** (ex.: AFP, GPC3, ALB, immune markers)
  antes e depois da correção.
- **Compare resultados com e sem ComBat** e reporte a sensibilidade.
- **Para scRNA-seq**, prefira métodos específicos de integração (Harmony,
  scVI, BBKNN, Seurat integration) em vez de ComBat gene a gene.
- **Para contagens**, use `ComBat-seq` ou modelos binomial negativo.
- **Evite ComBat quando lote e condição são perfeitamente confundidos.** Nesse
  caso, a conclusão honesta é que o efeito não é identificável; a solução é
  obter amostras adicionais ou usar comparações internas (ex.: tumor vs.
  adjacente dentro do TCGA).

---

## 2.6 Resumo em uma frase

O ComBat é uma ferramenta poderosa para remover efeitos técnicos quando o lote
é **separável** da biologia e quando as covariáveis biológicas são informadas;
na deconvolução celular, ele só deve ser usado com cautela, pois a composição
celular é simultaneamente biologia e possível fonte de lote — e corrigir o
lote sem preservar a composição pode apagar exatamente o sinal que se deseja
medir.

---

# 3. Como o `logFC` é calculado e interpretado no `limma`

> Documento de apoio ao seminário. Arquivo: `docs/03_logFC_limma.md`.
> Implementação: `scripts/02_expressao_diferencial.R`.
> Resultados: `resultados/tabelas/DE/`.

---

## 3.1 Ponto de partida: a escala logarítmica

No `limma`, o `logFC` (*log fold change*) é a **diferença de médias na escala
logarítmica**. Se a expressão já está em $\log_2(\text{TPM}+0{,}001)$, então:

$$
\text{logFC} = \overline{\log_2(\text{TPM}_{\text{grupo B}})}
             - \overline{\log_2(\text{TPM}_{\text{grupo A}})}
$$

Para voltar à escala de fold change:

$$
\text{FC} = 2^{\text{logFC}}
$$

Exemplos:

| logFC | FC aproximado | Interpretação |
|---:|---:|---|
| 0 | 1 | sem diferença |
| 0,585 | 1,5 | 50% maior no grupo B |
| 1 | 2 | 2× maior no grupo B |
| 2 | 4 | 4× maior no grupo B |
| −1 | 0,5 | metade no grupo B (2× menor) |
| −2 | 0,25 | 1/4 no grupo B |
| 6,97 | ≈125 | ≈125× maior (ex.: GPC3 no tumor) |
| −10,2 | ≈0,00085 | ≈1/1.177 (ex.: CLEC4M no tumor) |

**Importante:** o `logFC` do `limma` não é calculado a partir de contagens
brutas. O usuário deve fornecer uma matriz já em escala logarítmica (log2 TPM,
logCPM etc.) ou usar `voom` para converter contagens em logCPM com pesos. No
nosso caso, usamos log2(TPM + 0,001) do hub Toil/Xena; por isso o `limma` foi
aplicado diretamente, sem `voom`.

---

## 3.2 O modelo linear gene a gene

O `limma` ajusta, para **cada gene separadamente**, um modelo linear:

$$
\mathbf{y}_g = \mathbf{X}\boldsymbol{\beta}_g + \boldsymbol{\varepsilon}_g
$$

em que:

| Símbolo | Significado |
|---|---|
| $\mathbf{y}_g$ | vetor de expressão do gene $g$ nas $n$ amostras |
| $\mathbf{X}$ | matriz de desenho (*design matrix*), $n \times p$ |
| $\boldsymbol{\beta}_g$ | vetor de coeficientes do gene $g$ |
| $\boldsymbol{\varepsilon}_g$ | erros, assumidos com média zero e variância $\sigma_g^2$ |

A matriz $\mathbf{X}$ codifica o grupo biológico e as covariáveis. Cada coluna
é um coeficiente; cada linha é uma amostra.

### 3.2.1 Exemplo com dois grupos e covariáveis

Suponha o modelo:

```r
design <- model.matrix(~ 0 + grupo + sexo + idade, data = metadados)
```

com `grupo` tendo níveis `GTEx_Normal` e `TCGA_Tumor`. As colunas serão algo
como:

```text
(Intercept)  grupoTCGA_Tumor  sexomale  idade
```

Mas, com `~ 0 + grupo`, o `Intercept` é removido e as colunas ficam:

```text
grupoGTEx_Normal  grupoTCGA_Tumor  sexomale  idade
```

Nesse caso, `grupoGTEx_Normal` e `grupoTCGA_Tumor` são as **médias ajustadas**
de cada grupo (ajustadas por sexo e idade), e os coeficientes de sexo e idade
representam seus efeitos comuns aos dois grupos.

O `logFC` do contraste tumor − normal é:

$$
\widehat{\text{logFC}} =
\hat{\beta}_{\text{grupoTCGA\_Tumor}}
-
\hat{\beta}_{\text{grupoGTEx\_Normal}}
$$

Isso é feito pela matriz de contraste $\mathbf{c}$:

$$
\widehat{\text{logFC}} = \mathbf{c}^\top \hat{\boldsymbol{\beta}}_g
$$

No nosso código, o vetor `cvec` tem `-1` no coeficiente do grupo de referência
e `+1` no coeficiente do grupo de interesse. A função `contrasts.fit()` aplica
esse contraste a todos os genes de uma vez.

### 3.2.2 Equivalência com `makeContrasts`

O jeito mais explícito de fazer o mesmo é:

```r
design <- model.matrix(~ 0 + grupo + sexo + idade, data = metadados)
colnames(design) <- c("GTEx_Normal", "TCGA_Tumor", "sexo_male", "idade")

fit  <- lmFit(expr, design)
contr <- makeContrasts(
  TCGA_Tumor - GTEx_Normal,
  levels = design
)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

tt <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  sort.by = "P",
  adjust.method = "BH"
)
```

O resultado é o mesmo: um `logFC` por gene, já ajustado pelas covariáveis.

---

## 3.3 O que muda quando há covariáveis?

Sem covariáveis, o `logFC` é simplesmente a diferença das médias amostrais
entre os dois grupos. Com covariáveis, ele é a diferença das médias
**ajustadas** pelo modelo. Isso é importante quando:

- os grupos têm distribuições de sexo ou idade diferentes;
- há amostras desbalanceadas;
- existe um fator de confundimento que também afeta a expressão.

No nosso caso, a comparação TCGA tumor × GTEx normal ajusta por **sexo** e
**idade**, porque:

- a proporção de homens é alta nos dois grupos, mas não idêntica;
- a idade média difere (TCGA tumor ≈ 59,4 anos; GTEx ≈ 53,7 anos; Wilcoxon
  p = 3,96 × 10⁻⁶);
- sexo e idade afetam a expressão hepática e poderiam gerar diferenças
  espúrias se não fossem considerados.

Ainda assim, o ajuste por covariáveis **não corrige** o confundimento
estrutural entre estudo (TCGA vs. GTEx) e condição (tumor vs. normal). Ele
apenas remove parte da variação explicada por sexo e idade.

---

## 3.4 Análise pareada

No TCGA, 50 pacientes têm tumor e tecido adjacente. Nesse caso, podemos usar
o paciente como **fator de bloqueio**:

```r
design <- model.matrix(~ 0 + grupo + paciente)
```

Cada paciente tem duas amostras; o modelo estima um efeito fixo por paciente e
o contraste tumor − adjacente é calculado **dentro de cada paciente**. Isso
remove a variabilidade entre pacientes (idade, sexo, genética, etc.) e
geralmente aumenta o poder estatístico.

No nosso código, a análise pareada é a comparação
`TCGA_Tumor_vs_Adjacent_PAREADO`. Ela produziu 2.934 genes up e 776 down
(vs. 3.322 up e 808 down na análise não pareada), confirmando que a maior parte
do sinal é robusta.

---

## 3.5 `lmFit`, `eBayes` e `topTable`

### 3.5.1 `lmFit`

`lmFit(expr, design)` ajusta o modelo linear para cada gene. Ele retorna:

- coeficientes estimados;
- erros-padrão;
- variâncias residuais;
- graus de liberdade.

### 3.5.2 `eBayes`

`eBayes()` implementa o **moderated t-statistic** do limma. Em vez de estimar a
variância de cada gene isoladamente (o que é instável com poucas amostras), o
limma combina a variância do gene com uma variância *a priori* estimada a
partir de todos os genes. Isso “modera” os testes e aumenta o poder.

- `trend = TRUE`: a variância *a priori* pode depender da expressão média
  (genes de baixa expressão costumam ter maior variância).
- `robust = TRUE`: reduz a influência de genes outliers na estimativa da
  variância *a priori*.

### 3.5.3 `topTable`

`topTable()` monta a tabela final. As colunas são:

| Coluna | Significado |
|---|---|
| `logFC` | log2 fold change do contraste (grupo B − grupo A, conforme o contraste) |
| `AveExpr` | expressão média do gene no conjunto de amostras analisado (log2) |
| `t` | estatística t moderada do contraste |
| `P.Value` | p-valor nominal (não ajustado) |
| `adj.P.Val` | p-valor ajustado por múltiplas comparações (BH, no nosso caso) |
| `B` | log-odds de o gene estar diferencialmente expresso (evidência Bayesiana) |

Argumentos usados:

```r
topTable(
  fit2,
  coef      = 1,          # qual contraste usar
  number    = Inf,        # retornar todos os genes
  sort.by   = "P",        # ordenar pelo p-valor
  adjust.method = "BH"    # Benjamini-Hochberg
)
```

No nosso código, também adicionamos as colunas `gene`, `comparison`,
`direction` e `regulation` (Up/Down/NS).

---

## 3.6 Interpretação correta do sinal

O sinal do `logFC` depende de como o contraste foi definido. No nosso script:

```r
levels = c("GTEx_Normal", "TCGA_Tumor")
cvec["groupGTEx_Normal"] <- -1
cvec["groupTCGA_Tumor"] <-  1
```

Portanto:

- `logFC > 0` → expressão maior em **TCGA tumor** do que em **GTEx normal**;
- `logFC < 0` → expressão maior em **GTEx normal** (menor no tumor).

Para a comparação tumor vs. adjacente:

```r
levels = c("Adjacent", "Tumor")
cvec["groupAdjacent"] <- -1
cvec["groupTumor"]    <-  1
```

- `logFC > 0` → maior no tumor;
- `logFC < 0` → maior no tecido adjacente.

Sempre confira a coluna `direction` ou o subtítulo do volcano plot para não
inverter a interpretação.

---

## 3.7 Exemplos reais deste projeto

### Exemplo 1 — GPC3 (tumor TCGA vs. GTEx normal)

```text
gene   logFC   adj.P.Val   regulation
GPC3   6,969   1,24e-70    Up
```

- $2^{6,969} \approx 125$.
- GPC3 é ≈125× mais expresso no tumor do que no fígado normal.
- É um marcador clássico de hepatocarcinoma, o que dá validade biológica ao
  pipeline.

### Exemplo 2 — AFP

```text
AFP    2,852   3,23e-13    Up
```

- $2^{2,852} \approx 7,2$.
- AFP é ≈7,2× mais expresso no tumor.
- A magnitude é menor que a de GPC3, mas o gene é altamente significativo.

### Exemplo 3 — CLEC4M

```text
CLEC4M  -10,201  7,73e-93   Down
```

- $2^{-10,201} \approx 8,5 \times 10^{-4}$.
- Expressão ≈1.177× menor no tumor.
- CLEC4M é um receptor de células endoteliais sinusoidais hepáticas; sua perda
  é consistente com desdiferenciação/arquitetura vascular alterada no HCC.

### Exemplo 4 — XIST (GTEx masculino vs. feminino)

```text
XIST   -9,273   7,53e-50   Down
```

- $2^{-9,273} \approx 1,6 \times 10^{-3}$.
- XIST é ≈617× menor em homens.
- Esse é um **controle positivo** para a comparação de sexo: o pipeline
  recupera corretamente o gene de inativação do cromossomo X.

---

## 3.8 Erros comuns e boas práticas

1. **Calcular `logFC` a partir de TPM não transformado.** O `limma` espera
   escala log; se você passar TPM bruto, o `logFC` será uma diferença de TPM,
   não um log-razão. Use log2(TPM + pequeno offset) ou `voom`.
2. **Confundir `logFC` com fold change.** `logFC = 1` é 2×, não 1×.
3. **Ignorar o sinal do contraste.** `logFC` positivo/negativo depende de qual
   grupo está no numerador.
4. **Usar p-valor nominal em vez de FDR.** Com 12.401 genes, milhares de testes
   são feitos; use `adj.P.Val` (BH) e reporte o critério.
5. **Filtrar genes após ver o resultado.** O filtro de baixa expressão deve ser
   definido antes e aplicado igualmente a todas as comparações.
6. **Não ajustar por covariáveis importantes.** Sexo, idade e outros fatores
   podem gerar diferenças espúrias.
7. **Interpretar causalidade a partir de um `logFC`.** Expressão diferencial é
   associação, não prova de causalidade.
8. **Esquecer o confundimento de lote.** Um `logFC` grande entre TCGA e GTEx
   pode refletir biologia, lote ou ambos.

---

## 3.9 Resumo

- O `logFC` do `limma` é a diferença de médias na escala log2, obtida a partir
  de um modelo linear por gene.
- A matriz de desenho codifica grupos e covariáveis; o contraste define qual
  diferença será estimada.
- `eBayes` modera a variância e aumenta o poder; `topTable` resume os
  resultados.
- `2^logFC` converte o `logFC` em fold change.
- A interpretação correta exige conhecer o contraste, a escala dos dados e as
  covariáveis incluídas.

---

# 4. Volcano plots: construção e interpretação

> Arquivo: `docs/04_volcano_plots.md`.
> Script: `scripts/03_volcano_plots.R`.
> Figuras: `resultados/figuras/Volcano_*.png` (somente PNG).
> Tabelas: `resultados/tabelas/DE/` e `resultados/tabelas/ROTULOS_Volcano_*.csv`.

---

## 4.1 O que é um volcano plot?

O **volcano plot** combina duas informações de cada gene:

- **Eixo x:** `logFC` — magnitude e direção da diferença de expressão.
- **Eixo y:** `−log10(FDR)` — significância estatística ajustada para
  múltiplos testes. Quanto mais alto, menor o FDR.

Cada ponto é um gene. A forma de “vulcão” aparece porque genes com `logFC`
grande tendem a ter p-valores pequenos, mas genes com `logFC` pequeno podem ser
muito significativos se a variância for baixa.

### Regiões do gráfico

| Região | Interpretação |
|---|---|
| Direita superior | up-regulado e significativo |
| Esquerda superior | down-regulado e significativo |
| Centro (baixo) | sem evidência de diferença |
| Laterais inferiores | `logFC` grande, mas não significativo (alta variância/poucas amostras) |

### Critérios usados

- **FDR < 0,05** (Benjamini-Hochberg).
- **|logFC| > 1** (equivalente a fold change > 2 ou < 0,5).
- Cores:
  - vermelho = up-regulado;
  - azul = down-regulado;
  - cinza = não significativo.
- Linhas tracejadas: FDR = 0,05 e logFC = ±1.

> O critério |logFC| > 1 é uma escolha pragmática para priorizar efeitos
> biologicamente relevantes. Ele não é obrigatório; genes com |logFC| < 1
> podem ser significativos e devem ser consultados nas tabelas completas.

---

## 4.2 Os três volcano plots principais

### 4.2.1 LIHC (TCGA): tumor primário vs. tecido adjacente

![Volcano LIHC](resultados/figuras/Volcano_01_LIHC_TCGA_tumor_vs_adjacente.png)

*Figura V1. Volcano plot da comparação interna no TCGA: 369 tumores primários
vs. 50 tecidos adjacentes. Modelo `limma` com sexo e idade como covariáveis.*

**Resultado:** **3.322 genes up** e **808 genes down** no tumor
(FDR < 0,05 e |logFC| > 1).

**Genes up de destaque:** GABRD, PLVAP, DIPK2B, TOMM40L, HIGD1B, SLC26A6,
MSTO1, UBE2T, MXD3, CDC25C, BIRC5, GPC3, MMP9.

**Genes down de destaque:** CLEC4G, CLEC4M, FCN3, GDF2, CRHBP, CLEC1B, CFP,
MARCO, FCN2, ADAMTS13, CYP1A2, MT1H, SFRP5.

**Interpretação biológica:**

- **Perda de identidade hepática:** CLEC4G, CLEC4M, CLEC1B, FCN2, FCN3, GDF2,
  CRHBP, CYP1A2 e MARCO são genes de hepatócitos, células endoteliais
  sinusoidais ou células de Kupffer. A queda desses genes é consistente com
  desdiferenciação, perda de arquitetura lobular e alteração do
  microambiente imune.
- **Proliferação e ciclo celular:** UBE2T, CDC25C, BIRC5 (survivina), MXD3 e
  PLVAP indicam aumento de proliferação e remodelamento vascular.
- **Marcadores de HCC:** GPC3 e MMP9 sobem no tumor, validando a comparação.
- **GDF2 (BMP9):** reduzido no tumor; é um fator produzido por células
  estreladas hepáticas e endotélio, com papel na homeostase vascular.
- **SFRP5:** inibidor de Wnt, reduzido no tumor, compatível com ativação da via
  Wnt no HCC.

**Por que essa é a comparação interna mais valiosa?** Tumor e adjacente vêm do
mesmo estudo (TCGA), do mesmo paciente em 50 casos, e foram processados no
mesmo pipeline. Isso reduz o confundimento de lote, embora o adjacente possa
já apresentar alterações de campo de cancerização.

### 4.2.2 GTEx: masculino vs. feminino

![Volcano GTEx](resultados/figuras/Volcano_02_GTEx_masculino_vs_feminino.png)

*Figura V2. Volcano plot da comparação interna no GTEx (110 fígados normais):
77 homens vs. 33 mulheres. Modelo `limma` com idade como covariável.*

**Resultado:** **8 genes up** e **5 genes down** no sexo masculino
(FDR < 0,05 e |logFC| > 1).

**Genes up em homens:** DDX3Y, TTTY14, RPS4Y1, KDM5D, UTY, USP9Y, EIF1AY,
TSPAN8.

**Genes down em homens:** XIST, PZP, VIL1.

**Interpretação:**

- Os genes up são quase todos do **cromossomo Y** (DDX3Y, RPS4Y1, KDM5D, UTY,
  USP9Y, EIF1AY, TTTY14). Isso é um **controle positivo** de sexo: o pipeline
  recupera corretamente a diferença cromossômica.
- **XIST** é o gene da inativação do cromossomo X, expresso em mulheres e
  silenciado em homens; seu `logFC` negativo em homens é o esperado.
- O número pequeno de DEGs autossômicos reflete o poder limitado (n = 110) e o
  critério rigoroso. Genes como PZP e VIL1 podem ter relevância hepática, mas
  precisam de validação.
- **Por que essa comparação?** GTEx vs. GTEx não tem tumor vs. normal; a
  comparação interna mais informativa é por sexo (e por idade). Ela serve como
  controle de qualidade e mostra que o pipeline detecta sinal biológico real
  quando ele existe.

### 4.2.3 TCGA vs. GTEx: tumor vs. fígado normal

![Volcano TCGA vs GTEx](resultados/figuras/Volcano_03_TCGA_tumor_vs_GTEx_normal.png)

*Figura V3. Volcano plot da comparação entre estudos: 369 tumores TCGA vs.
110 fígados normais GTEx. Modelo `limma` com sexo e idade como covariáveis.*

**Resultado:** **2.245 genes up** e **1.143 genes down** no tumor
(FDR < 0,05 e |logFC| > 1).

**Genes up de destaque:** PLVAP, GBA1, RPL17P36, CD34, SEM1, PRPF6, GABRD,
LMNA, MDK, GPC3, BIRC5, VIM, MMP9, SPP1.

**Genes down de destaque:** PRSS53, CHRD, CLEC4M, ADAMTS13, RIPOR3, DNHD1,
NPIPB12, NPIPB5, GDF2, KRT13, CLEC1B, INS-IGF2, CYP1A2.

**Interpretação biológica:**

- **Ganho de proliferação, angiogênese e matriz:** PLVAP (endotélio
  fenestrado), CD34 (progenitores/endotélio), MDK (midquina), BIRC5, VIM,
  MMP9, SPP1.
- **Perda de identidade hepática:** CLEC4M, CLEC1B, GDF2, CYP1A2, CHRD,
  ADAMTS13.
- **GPC3 e AFP** sobem no tumor, reforçando a consistência com a literatura.
- **KRT13** e **INS-IGF2** aparecem com `logFC` muito negativo/positivo; alguns
  desses genes podem refletir diferenças de composição tecidual ou de
  processamento entre TCGA e GTEx.

**Cuidado crucial:** essa comparação é **confundida com o estudo**. Todo tumor
é TCGA; todo normal é GTEx. Portanto, o `logFC` estimado combina:

1. efeito biológico tumor vs. normal;
2. efeito de lote TCGA vs. GTEx;
3. diferenças de composição celular (tumor desmoplásico vs. fígado normal);
4. diferenças de seleção de doadores (idade, comorbidades, etc.).

A comparação **TCGA adjacente vs. GTEx normal** ajuda a estimar o componente
não tumoral: foram **584 up e 1.966 down** (20,6% dos genes), mostrando que os
“normais” dos dois estudos já diferem substancialmente. Isso não invalida a
comparação tumor vs. normal, mas exige moderação na interpretação dos efeitos
absolutos.

---

## 4.3 Volcano plots complementares

### Estágio III–IV vs. I–II (TCGA tumor)

![Volcano estágio](resultados/figuras/Volcano_05_TCGA_estagio_III_IV_vs_I_II.png)

**Resultado:** 3 up e 43 down (FDR < 0,05 e |logFC| > 1).

- **Up em III–IV:** MIR210HG, SLC6A8, CA9 (hipóxia).
- **Down em III–IV:** GNMT, CFHR4, CYP2C9, CFHR1, ACOT12, FTCD, HPX, SLC10A1.

**Interpretação:** tumores mais avançados perdem genes de metabolismo
hepático/diferenciação (GNMT, CYP2C9, ACOT12, FTCD, HPX, SLC10A1) e ganham
genes de hipóxia (CA9, MIR210HG). O número pequeno de DEGs reflete o tamanho
amostral e a heterogeneidade; ainda assim, o padrão é biologicamente coerente.

### Grau G3–G4 vs. G1–G2 (TCGA tumor)

![Volcano grau](resultados/figuras/Volcano_06_TCGA_grau_G3_G4_vs_G1_G2.png)

**Resultado:** 61 up e 111 down (FDR < 0,05 e |logFC| > 1).

- **Up em G3–G4:** HJURP, PTTG1, CCNB1, CDC20, TROAP, SKA3, AURKB, NUF2,
  KIF2C, CDKN3, KIFC1.
- **Down em G3–G4:** C6, CYP4F2, LINC01093, GHR, CCL14, CYP2A6, ESR1, GPD1,
  SLC10A1.

**Interpretação:** tumores de alto grau têm forte assinatura de **ciclo
celular/mitose** (AURKB, CCNB1, CDC20, NUF2, KIF2C) e perda de genes de
metabolismo/diferenciação. É um dos achados mais consistentes do ponto de vista
biológico.

### Sexo no TCGA (masculino vs. feminino)

![Volcano sexo TCGA](resultados/figuras/Volcano_07_TCGA_masculino_vs_feminino.png)

**Resultado:** 155 up e 64 down.

- Up em homens: genes do cromossomo Y (EIF1AY, DDX3Y, KDM5D, UTY, RPS4Y1,
  USP9Y, TTTY14).
- Down em homens: XIST e genes ligados ao X.

Esse resultado é um controle positivo e mostra que a comparação de sexo
funciona mesmo dentro do tumor.

### Idade no TCGA (≥65 vs. <65)

![Volcano idade TCGA](resultados/figuras/Volcano_08_TCGA_idade_65plus_vs_lt65.png)

**Resultado:** 43 up e 30 down. Efeito modesto; genes como LRRC37A7P, MFAP3L,
CES5A, SLC28A1, ALDH1L1 aparecem entre os mais significativos. A interpretação
é limitada pela natureza observacional e pela correlação entre idade e outras
variáveis clínicas.

### Idade no GTEx (≥60 vs. <60)

![Volcano idade GTEx](resultados/figuras/Volcano_09_GTEx_idade_60plus_vs_lt60.png)

**Resultado:** 2 up e 113 down. A maioria dos genes down em ≥60 está ligada a
ciclo celular/proliferação (KIF20A, TOP2A, RRM2) e endotélio (CD93, ECSCR).
Como o fígado adulto tem baixa proliferação, esse sinal pode refletir
diferenças de composição ou de viabilidade celular, e não necessariamente
envelhecimento hepático.

### TCGA adjacente vs. GTEx normal

![Volcano adjacente vs GTEx](resultados/figuras/Volcano_10_TCGA_adjacente_vs_GTEx_normal.png)

**Resultado:** 584 up e 1.966 down. É a comparação que mais evidencia o efeito
de estudo/lote em tecidos “normais”. Genes up no adjacente incluem APOC2, IGF2,
AMBP, CFP, FUCA1, ITIH2, CETP; genes down incluem CHRD, CCNL2, CLASRP,
SNRNP70, SLC17A9.

**Interpretação:** não se deve assumir que “adjacente TCGA” e “GTEx normal” são
intercambiáveis. O adjacente pode ter alterações de campo de cancerização,
inflamação e composição celular diferente; além disso, há diferenças de
processamento entre os estudos.

### Via KEGG hsa05417 (Lipídios e Aterosclerose)

![Volcano via LA](resultados/figuras/Volcano_11_LA_pathway_tumor_vs_GTEx.png)

**Resultado:** 14 up e 32 down entre os 213 genes da via.

- **Up no tumor:** BAX, HSPA4, HSP90AB1, PYCARD, PPARG, MMP1, LY96, ABCG1,
  MMP9, IKBKE, CD36, VCAM1, NOX1, PIK3R2.
- **Down no tumor:** CALML5, AGER, CXCL2, MIB2, CALML6, IKBKB.

**Interpretação:** a via de lipídios e aterosclerose está remodelada no HCC,
com ativação de componentes inflamatórios/oxidativos (PYCARD, LY96, VCAM1,
NOX1, MMP1/MMP9) e alteração de receptores de lipídios (CD36, ABCG1, PPARG).
A interpretação é limitada pelo pequeno número de genes da via e pelo
confundimento TCGA × GTEx.

### TCGA vs. GTEx após ComBat (demonstração auditada)

![Volcano ComBat](resultados/figuras/Volcano_12_TCGA_vs_GTEx_ComBat.png)

**Resultado:** a figura agora é um painel duplo: **antes do ComBat** (3.388
DEGs) e **depois do ComBat** (1 DEG, LINC01419). A versão original, com apenas
o painel “depois”, parecia vazia.

**Auditoria:** a figura original não estava tecnicamente em branco — tinha
≈3,65% de pixels não brancos (eixos, título e poucos pontos), mas 12.400 dos
12.401 genes ficaram amontoados perto da origem. A causa era científica: o
ComBat sem covariável removeu também a diferença tumor-normal, pois lote e
condição são perfeitamente confundidos. A versão corrigida tem ≈9,48% de
pixels não brancos e comunica explicitamente a sobrecorreção.

**Interpretação:** não é um resultado biológico; é uma demonstração de
**sobrecorreção**. Veja `docs/02_ComBat.md` e
`auditorias/AUDITORIA_VOLCANO_12.md` para a explicação completa.

---

## 4.4 Como os volcano plots foram gerados

Fluxo resumido:

1. Ajuste do modelo `limma` para cada comparação:
   ```r
   design <- model.matrix(~ 0 + grupo + sexo + idade, data = metadados)
   fit <- lmFit(expr, design)
   fit2 <- contrasts.fit(fit, contraste)
   fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
   tt <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
   ```
2. Classificação:
   ```r
   regulation <- ifelse(adj.P.Val < 0.05 & logFC >  1, "Up",
                 ifelse(adj.P.Val < 0.05 & logFC < -1, "Down", "NS"))
   ```
3. Volcano:
   ```r
   ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val), colour = regulation)) +
     geom_point() +
     geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
     geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
     geom_text_repel(data = rotulos, aes(label = gene))
   ```
4. Os genes rotulados são os DEGs com maior |logFC|, complementados por
   marcadores de interesse (AFP, GPC3, MMP9, VIM, CDH1, BIRC5 etc.).

### Por que usar FDR no eixo y?

Com 12.401 genes, usar o p-valor nominal geraria muitos falsos positivos. O
FDR (BH) controla a proporção esperada de falsos positivos entre os genes
chamados de significativos. No volcano, `−log10(FDR)` torna a significância
visualmente comparável entre comparações.

### Por que alguns genes conhecidos não aparecem rotulados?

O rótulo é limitado aos 10 DEGs com maior |logFC| + genes de interesse. Genes
com |logFC| menor, mas ainda significativos (por exemplo, AFP na comparação
tumor vs. adjacente), podem estar no gráfico sem rótulo. Consulte as tabelas
completas em `resultados/tabelas/DE/`.

---

## 4.5 Tabela-resumo dos contrastes

| Comparação | n testados | Up | Down | Total | % significativos |
|---|---:|---:|---:|---:|---:|
| TCGA tumor vs. TCGA adjacente | 12.401 | 3.322 | 808 | 4.130 | 33,3 |
| TCGA tumor vs. adjacente (pareado) | 12.401 | 2.934 | 776 | 3.710 | 29,9 |
| TCGA tumor vs. GTEx normal | 12.401 | 2.245 | 1.143 | 3.388 | 27,3 |
| TCGA adjacente vs. GTEx normal | 12.401 | 584 | 1.966 | 2.550 | 20,6 |
| TCGA masculino vs. feminino | 12.401 | 155 | 64 | 219 | 1,8 |
| TCGA G3–G4 vs. G1–G2 | 12.401 | 61 | 111 | 172 | 1,4 |
| GTEx ≥60 vs. <60 | 12.401 | 2 | 113 | 115 | 0,9 |
| TCGA ≥65 vs. <65 | 12.401 | 43 | 30 | 73 | 0,6 |
| TCGA estágio III–IV vs. I–II | 12.401 | 3 | 43 | 46 | 0,4 |
| Via LA: TCGA tumor vs. GTEx | 213 | 14 | 32 | 46 | 21,6 |
| Via LA: TCGA tumor vs. adjacente | 213 | 11 | 33 | 44 | 20,7 |
| GTEx masculino vs. feminino | 12.401 | 8 | 5 | 13 | 0,1 |

Critério: FDR < 0,05 e |logFC| > 1.

---

## 4.6 Limitações dos volcano plots

1. **Não mostram o tamanho amostral** nem a precisão de cada gene.
2. **O eixo y depende do ajuste de múltiplas comparações**; comparações com
   mais genes testados podem parecer “menos significativas”.
3. **Genes com baixa expressão** podem ter `logFC` instável; o filtro de
   expressão ajuda, mas não elimina o problema.
4. **Confundimento de lote** não aparece no gráfico: um ponto up pode ser
   biologia ou lote.
5. **O critério |logFC| > 1 é arbitrário**; para fins clínicos, podem ser
   relevantes efeitos menores.
6. **Volcano não substitui validação** por qPCR, IHC ou coortes independentes.

---

## 4.7 Resumo

- Os três volcano plots principais mostram:
  1. **LIHC TCGA** (tumor vs. adjacente): perda de identidade hepática e ganho
     de proliferação/remodelamento;
  2. **GTEx** (masculino vs. feminino): genes do cromossomo Y e XIST, um
     controle positivo;
  3. **TCGA vs. GTEx** (tumor vs. fígado normal): o maior conjunto de DEGs, mas
     confundido com estudo/lote.
- Os volcano complementares mostram que **grau** tem assinatura de ciclo
  celular, **estágio avançado** perde genes metabólicos e ganha hipóxia, e a
  **via LA** está remodelada no tumor.
- A interpretação biológica deve sempre considerar o desenho, as covariáveis e
  o confundimento.

---

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

**Outros grupos nos dados:** além de tumor/adjacente/GTEx, há 2 recidivas
(excluídas), 7 hepatocolangiocarcinomas mistos, 3 carcinomas fibrolamelares,
grupos virais (HBV/HCV/HBV+HCV) e 38 centros TCGA. Inventário em
`docs/10_grupos_alternativos.md`; análises separadas de cada grupo em
`docs/11_grupos_separados.md`.

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

---

# 6. Dicionário de dados e arquivos

> Arquivo: `docs/06_dicionario_dados.md`.

---

## 6.1 Fontes dos dados

| Arquivo | Fonte | URL |
|---|---|---|
| `dados/raw/LIHC_clinicalMatrix` | UCSC Xena / TCGA LIHC clinical matrix | `https://tcga.xenahubs.net/download/TCGA.LIHC.sampleMap/LIHC_clinicalMatrix` |
| `dados/raw/TCGA-LIHC.clinical.tsv.gz` | UCSC Xena / GDC TCGA-LIHC | `https://gdc.xenahubs.net/download/TCGA-LIHC.clinical.tsv.gz` |
| `dados/raw/GTEX_phenotype.gz` | UCSC Xena / GTEx phenotype | `https://toil.xenahubs.net/download/GTEX_phenotype.gz` |
| `dados/raw/TcgaTargetGTEX_phenotype.txt.gz` | UCSC Xena / Toil phenotype | `https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz` |
| `dados/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` | GTEx Portal v8 | `https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` |
| `dados/raw/liver.tsv` | painel hepático (213 genes da via KEGG hsa05417 + metadados); cópia local do projeto principal | repositório `hepatocellular-carcinoma` |
| `TcgaTargetGtex_rsem_gene_tpm.gz` | UCSC Xena / Toil (transcriptoma completo) | `https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz` |

O transcriptoma completo (≈1,3 GB) **não** foi versionado dentro de
`seminários_ryan/`; ele foi processado e a matriz filtrada (17 MB) foi salva em
`dados/processed/expressao_log2tpm_filtrada.tsv.gz`. Para reproduzir do zero,
use `scripts/00_download_dados.R` ou baixe manualmente pela URL acima e defina
a variável de ambiente `TCGA_GTEX_TPM` com o caminho do arquivo.

---

## 6.2 Seleção das amostras

- **TCGA-LIHC tumor primário:** `_sample_type == "Primary Tumor"` e
  `_study == "TCGA"` → 369 amostras.
- **TCGA-LIHC tecido adjacente:** `_sample_type == "Solid Tissue Normal"` →
  50 amostras.
- **GTEx fígado normal:** `_study == "GTEX"`, `_primary_site == "Liver"`,
  `_sample_type == "Normal Tissue"` → 110 amostras.
- **Recidivas:** 2 amostras `Recurrent Tumor` (`TCGA-DD-AACA-02` e
  `TCGA-ZS-A9CF-02`) foram **excluídas** da análise principal.
- **Total:** 531 amostras no metadado completo; 529 na análise principal.

---

## 6.3 `metadata_amostras.tsv` e `metadata_analise.tsv`

| Coluna | Descrição |
|---|---|
| `sample` | identificador da amostra (TCGA barcode ou GTEx sample ID) |
| `samples` | cópia do identificador (coluna original do liver.tsv) |
| `detailed_category` | categoria detalhada do Xena (`Liver Hepatocellular Carcinoma`, etc.) |
| `DFI`, `DFI.time` | doença livre de intervalo e tempo (quando disponível) |
| `DSS`, `DSS.time` | sobrevida específica da doença e tempo |
| `OS`, `OS.time` | sobrevida global e tempo (dias) |
| `PFI`, `PFI.time` | intervalo livre de progressão e tempo |
| `sex_liver` | sexo original do liver.tsv (`_gender`) |
| `sample_type` | tipo de amostra (`Primary Tumor`, `Solid Tissue Normal`, `Normal Tissue`, `Recurrent Tumor`) |
| `study` | `TCGA` ou `GTEX` |
| `main_category` | categoria principal do Xena (`TCGA Liver Hepatocellular Carcinoma` ou `GTEX Liver`) |
| `condition` | `Tumor`, `Adjacent`, `Normal` ou `Recurrent` |
| `group` | `TCGA_Tumor`, `TCGA_Adjacent`, `GTEx_Normal`, `TCGA_Recurrent` |
| `dataset` | `TCGA-LIHC` ou `GTEx-Liver` |
| `subject` | identificador do paciente/doador (TCGA: `TCGA-XX-YYYY`; GTEx: `GTEX-XXXX`) |
| `sex_tcga` | sexo do clínico TCGA |
| `age` | idade em anos (TCGA: idade ao diagnóstico; GTEx: ponto médio da faixa) |
| `stage_raw` | estágio AJCC original (texto) |
| `stage` | estágio harmonizado (`I`, `II`, `III`, `IV`) |
| `grade_raw` | grau histológico original (`G1`–`G4`, `Not Reported`) |
| `grade` | grau harmonizado (`G1`–`G4`) |
| `histology` | tipo histológico |
| `vital_status` | status vital |
| `days_to_death` | dias até o óbito |
| `days_to_last_followup` | dias de seguimento |
| `fibrosis_ishak` | escore de fibrose de Ishak |
| `afp` | alfa-fetoproteína (valor laboratorial) |
| `viral_hepatitis` | sorologia para hepatites virais |
| `child_pugh` | classificação de Child-Pugh |
| `sex_gtex` | sexo do GTEx (1 = masculino, 2 = feminino) |
| `age_gtex` | idade do GTEx (ponto médio da faixa) |
| `sex` | sexo harmonizado (`female`/`male`) |
| `age_group` | `<65` ou `>=65` |
| `age_group60` | `<60` ou `>=60` |
| `stage_group` | `I-II` ou `III-IV` |
| `grade_group` | `G1-G2` ou `G3-G4` |
| `tss` | código do centro TCGA (Tissue Source Site) |
| `viral_status` | `None/Other`, `HBV`, `HCV`, `HBV+HCV` ou `Unknown` |
| `afp_log2` | log2(AFP + 1) |
| `afp_high` | AFP `<400` ou `>=400` ng/mL |
| `afp_group` | AFP `<20`, `20-399` ou `>=400` |
| `fibrosis_group` | Ishak `0`, `1-4`, `5` ou `6` |
| `fibrosis_advanced` | Ishak `<5` ou `>=5` |
| `child_pugh_group` | Child-Pugh `A` ou `B-C` |
| `histology_group` | `HCC` ou `Other` |
| `vascular_invasion` | `None/Unknown`, `Micro` ou `Macro` |
| `race`, `race_group` | raça declarada e agrupada |
| `ethnicity`, `ethnicity_group` | etnia declarada e agrupada |
| `pathologic_t`, `pathologic_n`, `pathologic_m` | T/N/M patológicos |
| `ajcc_t`, `ajcc_n`, `ajcc_m` | T/N/M AJCC (GDC) |
| `residual_tumor` | R0, R1, R2 ou RX |
| `neoplasm_status` | status neoplásico no seguimento |

> Observação: para as 50 amostras adjacentes, `stage` e `grade` refletem o
> tumor do mesmo paciente, não o tecido adjacente em si. As análises de
> estágio/grau usam apenas `group == "TCGA_Tumor"`.

---

## 6.4 Matrizes de expressão

### `expressao_log2tpm_filtrada.tsv.gz`

- Formato: TSV gzipado, primeira coluna `gene` (símbolo HUGO), demais colunas
  amostras.
- Valores: log2(TPM + 0,001).
- Dimensões: 12.401 genes × 529 amostras.
- Filtro: expressão ≥ 1 em ≥ 20% das amostras; variância > 0.

### `expressao_LA_pathway.tsv.gz`

- Formato: TSV gzipado, primeira coluna `gene`.
- Genes: 213 genes da via KEGG hsa05417 (Lipídios e Aterosclerose).
- Valores: log2(TPM + 0,001).
- Amostras: 531 (inclui recidivas, mas as análises principais usam 529).

### `resumo_preparacao.txt`

- Resumo textual da etapa 00, com contagens, dimensões e sessão R.

---

## 6.5 Tabelas de expressão diferencial

Arquivos em `resultados/tabelas/DE/`:

| Arquivo | Conteúdo |
|---|---|
| `DEG_<comparacao>_completo.csv` | todos os genes testados, com estatísticas completas |
| `DEG_<comparacao>_significativos.csv` | apenas DEGs (FDR < 0,05 e \|logFC\| > 1) |
| `TOP20_<comparacao>.csv` | top 20 up + top 20 down |
| `RESUMO_DEGs_TODAS_COMPARACOES.csv` | contagens de up/down por comparação |
| `COMPARACAO_ComBat_antes_depois.csv` | comparação dos DEGs antes/depois do ComBat |

### Colunas das tabelas de DEG

| Coluna | Descrição |
|---|---|
| `gene` | símbolo HUGO |
| `comparison` | nome da comparação |
| `direction` | direção do contraste (ex.: `TCGA_Tumor > GTEx_Normal`) |
| `logFC` | log2 fold change |
| `AveExpr` | expressão média do gene (log2) |
| `t` | estatística t moderada do limma |
| `P.Value` | p-valor nominal |
| `adj.P.Val` | p-valor ajustado por BH (FDR) |
| `B` | log-odds bayesiano de diferença |
| `regulation` | `Up`, `Down` ou `NS` |

### Nomes das comparações

- `TCGA_Tumor_vs_TCGA_Adjacent`
- `TCGA_Tumor_vs_Adjacent_PAREADO`
- `TCGA_Tumor_vs_GTEx_Normal`
- `TCGA_Adjacent_vs_GTEx_Normal`
- `GTEx_Male_vs_Female`
- `GTEx_Age60plus_vs_lt60`
- `TCGA_Stage_III_IV_vs_I_II`
- `TCGA_Grade_G3_G4_vs_G1_G2`
- `TCGA_Male_vs_Female`
- `TCGA_Age65plus_vs_lt65`
- `LApathway_TCGA_Tumor_vs_GTEx_Normal`
- `LApathway_TCGA_Tumor_vs_TCGA_Adjacent`
- `TCGA_Tumor_vs_GTEx_Normal_ComBat` (demonstração de sobrecorreção)

---

## 6.6 Tabelas clínicas avançadas

| Arquivo | Conteúdo |
|---|---|
| `13_fenotipos_missingness.csv` | missingness das variáveis clínicas |
| `14_fenotipos_testes.csv` | testes de associação e Cramér's V |
| `15_AFP_por_subgrupo.csv` | resumo de AFP por estágio, grau, sexo e vírus |
| `16_KM_logrank.csv` | log-rank das curvas de sobrevida |
| `17_cox_univariado.csv` | Cox univariado |
| `18_cox_multivariado.csv` | Cox multivariado (Modelos A e B) |
| `19_correlacao_clinica.csv` | correlações de Spearman |
| `20_cox_ph_assumption.csv` | `cox.zph` (riscos proporcionais) |
| `21_inventario_grupos.csv` | inventário de tipos de amostra, histologia, vírus e TSS |
| `22_histologia_marcadores.csv` | medianas de marcadores por subtipo histológico |
| `23_recidivas_clinico.csv` | dados clínicos das 2 recidivas |
| `24_recidivas_comparacao.csv` | recidivas vs. resumo dos tumores primários |
| `25_recidivas_expressao.csv` | genes mais extremos das recidivas no painel LA |
| `26_recidivas_correlacao.csv` | similaridade global das recidivas com cada grupo |
| `27_histologia_clinico.csv` | resumo clínico por subtipo histológico |
| `28_histologia_testes.csv` | Kruskal-Wallis e comparações de marcadores por histologia |
| `29_viral_clinico.csv` | resumo clínico por grupo viral |
| `30_viral_testes.csv` | testes entre grupos virais e conhecido vs. Unknown |
| `31_viral_marcadores.csv` | medianas de marcadores por grupo viral |
| `32_tss_resumo.csv` | resumo clínico por centro TCGA |
| `33_tss_testes.csv` | testes entre centros com n ≥ 10 |
| `34_tss_marcadores.csv` | medianas de marcadores por TSS |
| `35_tss_umap.csv` | associação entre UMAP e TSS |

## 6.7 Figuras

| Arquivo | Conteúdo |
|---|---|
| `01_amostras_por_grupo.png` | distribuição das amostras |
| `02_sexo_por_grupo.png` | composição por sexo |
| `03_idade_por_grupo.png` | idade por grupo (violino + boxplot) |
| `04_estagio_TCGA.png` | estágio AJCC |
| `05_grau_TCGA.png` | grau histológico |
| `06_estagio_x_grau_TCGA.png` | estágio × grau |
| `07_UMAP_grupo.png` | UMAP por grupo |
| `08_UMAP_estudo.png` | UMAP por estudo |
| `09_UMAP_sexo_TCGA.png` | UMAP por sexo (TCGA tumor) |
| `10_UMAP_estagio_grau_TCGA.png` | UMAP por estágio/grau |
| `11_heatmap_top50_variaveis.png` | heatmap top 50 genes variáveis |
| `12_violino_marcadores.png` | violinos + boxplot de marcadores |
| `13_UMAP_ComBat_antes_depois.png` | UMAP antes/depois do ComBat |
| `14_AFP_ComBat_antes_depois.png` | AFP antes/depois do ComBat (violino + boxplot) |
| `15_AFP_por_estagio_grau.png` | AFP por estágio e grau (violino + boxplot) |
| `16_missingness_clinico.png` | missingness das variáveis clínicas |
| `17_KM_OS_estagio.png` | Kaplan-Meier por estágio |
| `18_KM_OS_grau.png` | Kaplan-Meier por grau |
| `19_KM_OS_sexo.png` | Kaplan-Meier por sexo |
| `20_KM_OS_idade.png` | Kaplan-Meier por idade |
| `21_KM_OS_AFP.png` | Kaplan-Meier por AFP |
| `22_KM_OS_fibrose.png` | Kaplan-Meier por fibrose |
| `23_KM_OS_painel.png` | painel de curvas KM |
| `24_forest_cox_multivariado.png` | forest plot do Cox multivariado |
| `25_correlacao_clinica.png` | heatmap de correlações clínicas |
| `26_grupos_dados.png` | tipos de amostra e subtipos histológicos |
| `27_histologia_marcadores.png` | marcadores por subtipo histológico |
| `28_recidivas_expressao.png` | genes extremos das recidivas no painel LA |
| `29_recidivas_correlacao.png` | correlação global das recidivas com os grupos |
| `30_histologia_clinico.png` | idade e AFP por subtipo histológico |
| `31_viral_clinico.png` | AFP por grupo viral |
| `32_viral_marcadores.png` | marcadores por grupo viral |
| `33_tss_distribuicao.png` | amostras por centro TCGA |
| `34_tss_umap.png` | UMAP colorido por TSS |
| `35_tss_marcadores_heatmap.png` | medianas de marcadores por TSS |
| `Volcano_01_LIHC_TCGA_tumor_vs_adjacente.png` | volcano principal LIHC |
| `Volcano_02_GTEx_masculino_vs_feminino.png` | volcano principal GTEx |
| `Volcano_03_TCGA_tumor_vs_GTEx_normal.png` | volcano principal TCGA vs. GTEx |
| `Volcano_04...12...` | volcanos complementares |

---

## 6.8 Convenções

- **logFC:** sempre em log2; positivo = maior no grupo do numerador do
  contraste.
- **FDR:** Benjamini-Hochberg.
- **Significância:** FDR < 0,05 e |logFC| > 1, salvo indicação contrária.
- **Idade:** anos; no GTEx, ponto médio da faixa.
- **Estágio:** AJCC harmonizado em I–IV.
- **Grau:** G1–G4.
- **Cores dos grupos:** TCGA tumor (vermelho), TCGA adjacente (laranja),
  GTEx normal (azul).

---

# 7. Guia de reprodutibilidade

> Arquivo: `docs/07_reproducibilidade.md`.

---

## 7.1 Requisitos

- **R ≥ 4.5** (testado com R 4.6.1 no Windows 11).
- Pacotes listados em `requirements_R.txt`.
- Espaço em disco:
  - ~20 MB para a pasta `seminários_ryan/` (matriz processada incluída);
  - ~1,3 GB adicionais se baixar o transcriptoma completo do Toil.
- Memória: ≥ 8 GB recomendado para ler o TPM completo; a matriz filtrada
  (12.401 × 529) usa bem menos.

Verifique a sessão usada:

```bash
cat docs/sessionInfo.txt
```

---

## 7.2 Ordem de execução

```text
00_download_dados.R        (opcional; baixa dados brutos)
        │
        ▼
00_preparar_dados.R        (metadados + matriz log2 TPM)
        │
        ▼
01_analise_exploratoria.R  (fenótipos, UMAP, heatmap, marcadores)
        │
        ▼
02_expressao_diferencial.R (limma em todas as comparações)
        │
        ▼
03_volcano_plots.R         (12 volcano plots)
        │
        ▼
04_ComBat_demonstracao.R   (sensibilidade/limitação do ComBat)
        │
        ▼
05_fenotipos_avancados.R   (missingness, testes, KM, Cox)
        │
        ▼
06_grupos_alternativos.R   (recidivas, histologia, vírus, TSS)
        │
        ▼
07_grupos_separados.R      (análises separadas de cada grupo)
        │
        ▼
10_auditoria_dados_e_volcano.R  (auditorias)
11_benchmarks.R                 (benchmarks)
12_compliance.R                 (checklist final)
        │
        ▼
09_gerar_relatorio_md.py        (opcional; Markdown consolidado)
```

`run_all.R` executa essa sequência automaticamente.

---

## 7.3 Passo a passo

### 7.3.1 Entrar na pasta

```bash
cd "caminho/para/hepatocellular-carcinoma/seminários_ryan"
```

### 7.3.2 (Opcional) Baixar dados

Os arquivos pequenos já estão em `dados/raw/`. Para rebaixar:

```bash
Rscript scripts/00_download_dados.R
```

Para incluir o transcriptoma completo (~1,3 GB):

```bash
BAIXAR_TRANSCRIPTOMA_COMPLETO=TRUE Rscript scripts/00_download_dados.R
```

O script baixa:

- `LIHC_clinicalMatrix` (TCGA, Xena);
- `TCGA-LIHC.clinical.tsv.gz` (GDC, Xena);
- `GTEX_phenotype.gz` e `TcgaTargetGTEX_phenotype.txt.gz` (Toil/Xena);
- `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` (GTEx Portal);
- opcionalmente `TcgaTargetGtex_rsem_gene_tpm.gz`.

### 7.3.3 Preparar dados

Se você baixou o TPM completo:

```bash
TCGA_GTEX_TPM="C:/caminho/TcgaTargetGtex_rsem_gene_tpm.gz" \
  Rscript scripts/00_preparar_dados.R
```

Se o TPM estiver em `$TEMP/hepato_full/TcgaTargetGtex_rsem_gene_tpm.gz`
(caminho usado no desenvolvimento), o script encontra automaticamente.

Se o TPM completo não estiver disponível, o script usa o painel
`../dados/raw/liver.tsv` (213 genes da via LA) e gera
`expressao_LA_pathway.tsv.gz`. As etapas seguintes usam a matriz completa se
ela já existir em `dados/processed/expressao_log2tpm_filtrada.tsv.gz`; caso
contrário, usam o painel LA.

**Saídas esperadas:**

- `dados/processed/metadata_amostras.tsv`
- `dados/processed/metadata_analise.tsv`
- `dados/processed/expressao_log2tpm_filtrada.tsv.gz` (se TPM completo)
- `dados/processed/expressao_LA_pathway.tsv.gz`
- `dados/processed/resumo_preparacao.txt`

### 7.3.4 Análise exploratória

```bash
Rscript scripts/01_analise_exploratoria.R
```

**Saídas:** `resultados/tabelas/01–12_*.csv`,
`resultados/figuras/01–12_*.png` e `resultados/qc/resumo_analise_exploratoria.txt`.

### 7.3.5 Expressão diferencial

```bash
Rscript scripts/02_expressao_diferencial.R
```

**Saídas:** `resultados/tabelas/DE/DEG_*_completo.csv`,
`DEG_*_significativos.csv`, `TOP20_*`, `RESUMO_DEGs_TODAS_COMPARACOES.csv`,
`resultados/qc/resumo_expressao_diferencial.txt`.

> Observação: nessa etapa, o ComBat com covariável biológica falha
> propositalmente (confundimento TCGA × GTEx). A mensagem é capturada e
> registrada; o pipeline continua. A demonstração completa está no script 04.

### 7.3.6 Volcano plots

```bash
Rscript scripts/03_volcano_plots.R
```

**Saídas:** `resultados/figuras/Volcano_01...12...png` (somente PNG) e
`resultados/tabelas/ROTULOS_Volcano_*.csv`.

### 7.3.7 ComBat

```bash
Rscript scripts/04_ComBat_demonstracao.R
```

**Saídas:**

- `resultados/tabelas/DE/DEG_TCGA_Tumor_vs_GTEx_Normal_ComBat_completo.csv`
  (resultado após ComBat sem covariável);
- `resultados/tabelas/DE/COMPARACAO_ComBat_antes_depois.csv`;
- `resultados/figuras/13_UMAP_ComBat_antes_depois.png`;
- `resultados/figuras/14_AFP_ComBat_antes_depois.png`;
- `resultados/qc/resumo_ComBat_demonstracao.txt`.

> Esse script sobrescreve `DEG_TCGA_Tumor_vs_GTEx_Normal_completo.csv` com a
> mesma análise limma (antes do ComBat), garantindo consistência com o script
> 02.

### 7.3.8 Fenótipos avançados

```bash
Rscript scripts/05_fenotipos_avancados.R
```

**Saídas:** `resultados/tabelas/13–20_*.csv`,
`resultados/figuras/15–25_*.png` e
`resultados/qc/resumo_fenotipos_avancados.txt`.

### 7.3.9 Grupos alternativos

```bash
Rscript scripts/06_grupos_alternativos.R
```

**Saídas:** `resultados/tabelas/21–22_*.csv`,
`resultados/figuras/26–27_*.png` e
`resultados/qc/resumo_grupos_alternativos.txt`.

### 7.3.10 Grupos separados

```bash
Rscript scripts/07_grupos_separados.R
```

**Saídas:** `resultados/tabelas/23–35_*.csv`,
`resultados/figuras/28–35_*.png` e
`resultados/qc/resumo_grupos_separados.txt`.

### 7.3.11 Auditorias, benchmarks e compliance

```bash
Rscript scripts/10_auditoria_dados_e_volcano.R
Rscript scripts/11_benchmarks.R
Rscript scripts/12_compliance.R
```

**Saídas:** `auditorias/AUDITORIA_DADOS.md`,
`auditorias/AUDITORIA_VOLCANO_12.md`, `auditorias/BENCHMARKS.md`,
`auditorias/COMPLIANCE.md` e respectivos CSVs.

### 7.3.12 Relatório Markdown consolidado (opcional)

```bash
python scripts/09_gerar_relatorio_md.py
```

---

## 7.4 Execução completa

```bash
Rscript scripts/run_all.R
```

O `run_all.R`:

1. detecta o TPM completo em `TCGA_GTEX_TPM` ou em
   `$TEMP/hepato_full/TcgaTargetGtex_rsem_gene_tpm.gz`;
2. roda 00 → 01 → 02 → 03;
3. roda 04 somente se a matriz completa existir;
4. roda 05 (fenótipos avançados);
5. roda 06 (grupos alternativos);
6. roda 07 (análises separadas dos grupos alternativos);
7. roda 10, 11 e 12 (auditorias, benchmarks e compliance);
8. para em caso de erro.

### 7.4.1 (Opcional) Gerar Markdown consolidado

```bash
python scripts/09_gerar_relatorio_md.py
```

Não requer dependências externas. Gera `relatorio_completo.md` na raiz,
concatenando os documentos de `docs/` e ajustando os caminhos das figuras.
O relatório é **somente Markdown**; não há versão HTML.

---

## 7.5 Estrutura de saídas

```text
resultados/
├── figuras/
│   ├── 01_amostras_por_grupo.png
│   ├── 02_sexo_por_grupo.png
│   ├── 03_idade_por_grupo.png
│   ├── 04_estagio_TCGA.png
│   ├── 05_grau_TCGA.png
│   ├── 06_estagio_x_grau_TCGA.png
│   ├── 07_UMAP_grupo.png
│   ├── 08_UMAP_estudo.png
│   ├── 09_UMAP_sexo_TCGA.png
│   ├── 10_UMAP_estagio_grau_TCGA.png
│   ├── 11_heatmap_top50_variaveis.png
│   ├── 12_violino_marcadores.png
│   ├── 13_UMAP_ComBat_antes_depois.png
│   ├── 14_AFP_ComBat_antes_depois.png
│   ├── 15_AFP_por_estagio_grau.png
│   ├── 16_missingness_clinico.png
│   ├── 17–23_KM_OS_*.png
│   ├── 24_forest_cox_multivariado.png
│   ├── 25_correlacao_clinica.png
│   ├── 26_grupos_dados.png
│   ├── 27_histologia_marcadores.png
│   ├── 28_recidivas_expressao.png
│   ├── 29_recidivas_correlacao.png
│   ├── 30_histologia_clinico.png
│   ├── 31_viral_clinico.png
│   ├── 32_viral_marcadores.png
│   ├── 33_tss_distribuicao.png
│   ├── 34_tss_umap.png
│   ├── 35_tss_marcadores_heatmap.png
│   └── Volcano_01...12...png
├── tabelas/
│   ├── 01–35_*.csv
│   └── DE/
│       ├── DEG_*_completo.csv
│       ├── DEG_*_significativos.csv
│       ├── TOP20_*.csv
│       ├── RESUMO_DEGs_TODAS_COMPARACOES.csv
│       └── COMPARACAO_ComBat_antes_depois.csv
└── qc/
    ├── log_00_preparar_dados.txt
    ├── log_01_analise_exploratoria.txt
    ├── log_02_expressao_diferencial.txt
    ├── log_03_volcano.txt
    ├── log_04_ComBat.txt
    ├── log_05_fenotipos_avancados.txt
    ├── log_06_grupos_alternativos.txt
    ├── resumo_grupos_separados.txt
    ├── log_10_auditoria.txt
    ├── log_11_benchmarks.txt
    ├── log_12_compliance.txt
    ├── resumo_analise_exploratoria.txt
    ├── resumo_expressao_diferencial.txt
    ├── resumo_ComBat_demonstracao.txt
    ├── UMAP_ComBat_coordenadas.csv
    └── resumo_fenotipos_avancados.txt
```

Os relatórios de auditoria, benchmark e compliance ficam em `auditorias/`:

```text
auditorias/
├── AUDITORIA_DADOS.md
├── AUDITORIA_VOLCANO_12.md
├── BENCHMARKS.md
├── COMPLIANCE.md
├── auditoria_dados_resumo.csv
├── auditoria_volcano12_diagnostico.csv
├── BENCHMARKS_tabela.csv
└── COMPLIANCE_tabela.csv
```

---

## 7.6 Solução de problemas

### “O arquivo TPM não é encontrado”

- Verifique se `TCGA_GTEX_TPM` aponta para o arquivo `.gz`.
- No Windows, use barras normais (`C:/...`) ou barras invertidas escapadas.
- Se não tiver o TPM completo, rode sem ele: o script 00 usa o painel LA. Se a
  matriz completa já existir em `dados/processed/expressao_log2tpm_filtrada.tsv.gz`,
  as etapas seguintes continuam usando essa matriz.

### “The covariate is confounded with batch!”

- É o comportamento **esperado** no script 04 para `mod = ~ grupo`.
- O script captura o erro e segue para a demonstração sem covariável.
- Não é um bug; é a limitação central do ComBat neste desenho.

### “Poucas amostras” ou “NA”

- Confira `metadata_analise.tsv`; amostras sem idade/sexo podem ser removidas
  por `complete.cases` nas comparações ajustadas.
- O TCGA tem 1 tumor sem idade; o GTEx tem idade para todos os 110.

### Acentuação em caminhos

- A pasta `seminários_ryan` tem acento. Em sistemas com locale não UTF-8,
  renomeie para `seminarios_ryan` ou execute a partir da pasta com
  `cd`. Os scripts usam caminhos relativos e `normalizePath`, então funcionam
  no Windows 11 com locale pt-BR.

### Diferenças numéricas

- Pequenas diferenças podem ocorrer entre versões de R/Bioconductor.
- As versões usadas estão em `docs/sessionInfo.txt`.
- O número de genes após o filtro pode variar se o mapeamento
  Ensembl → símbolo mudar.

---

## 7.7 Integridade e versionamento

- A pasta contém apenas resultados leves (figuras, tabelas, matriz processada).
- O TPM bruto de 1,3 GB **não** deve ser versionado.
- Se for usar Git, adicione ao `.gitignore`:

```gitignore
seminários_ryan/dados/raw/TcgaTargetGtex_rsem_gene_tpm.gz
seminários_ryan/resultados/qc/*.log
```

- Os scripts são determinísticos, exceto:
  - `ComBat` (usa `set.seed(2026)` no script 04);
  - `ggrepel` (posicionamento de rótulos; `seed = 2026` no volcano);
  - pequenas variações de versão de pacote.

---

# 8. Fenótipos clínicos avançados, sobrevida e Cox

> Arquivo: `docs/08_fenotipos_avancados.md`.
> Script: `scripts/05_fenotipos_avancados.R`.
> Saídas: `resultados/tabelas/13–20_*.csv`, `resultados/figuras/15–25_*.png`.

---

## 8.1 Objetivo

A análise exploratória inicial descreveu sexo, idade, estágio e grau. Esta
etapa amplia a exploração clínica do TCGA-LIHC para:

- quantificar **missingness** e qualidade das variáveis;
- testar associações entre fenótipos;
- explorar **AFP**, **fibrose**, **sorologia viral**, **Child-Pugh**,
  **histologia**, **invasão vascular** e **raça/etnia**;
- analisar **sobrevida global** com Kaplan-Meier e log-rank;
- ajustar modelos de **Cox univariados e multivariados**;
- verificar a suposição de **riscos proporcionais** (`cox.zph`);
- medir correlações entre variáveis clínicas.

Todos os resultados são **exploratórios** e observacionais.

---

## 8.2 Missingness

![Missingness clínico](resultados/figuras/16_missingness_clinico.png)

*Figura 16. Percentual de valores ausentes/desconhecidos nas variáveis clínicas
dos 369 tumores primários TCGA.*

| Variável | Ausente/desconhecido | % |
|---|---:|---:|
| Invasão vascular | 261 | 70,7 |
| Sorologia viral | 206 | 55,8 |
| Fibrose (Ishak) | 159 | 43,1 |
| Fibrose avançada | 159 | 43,1 |
| Child-Pugh | 132 | 35,8 |
| AFP | 93 | 25,2 |
| AFP ≥400 | 93 | 25,2 |
| Estágio | 24 | 6,5 |
| Etnia | 19 | 5,1 |
| Raça | 10 | 2,7 |
| Grau | 5 | 1,4 |
| Idade | 1 | 0,3 |
| OS/time | 1 | 0,3 |
| Sexo, histologia | 0 | 0 |

**Implicação:** as análises de AFP, fibrose, Child-Pugh, vírus e invasão
vascular têm poder limitado e podem sofrer viés de seleção (amostras com
informação disponível podem não ser representativas).

---

## 8.3 Associações entre fenótipos

A tabela completa está em
`resultados/tabelas/14_fenotipos_testes.csv` (24 testes, com p nominal, p
ajustado por BH e V de Cramér). Para as variáveis **vírus** e **raça**, a
categoria `Unknown` foi excluída dos testes, pois não representa um grupo
biológico.

### Associações significativas (BH < 0,05)

| Associação | n | Teste | p | V de Cramér |
|---|---:|---|---:|---:|
| Child-Pugh × estágio | 345 | Qui-quadrado | 5,5 × 10⁻⁸ | 0,293 |
| log2(AFP+1) × grau | 275 | Wilcoxon | 4,3 × 10⁻⁶ | — |
| AFP ≥400 × grau | 275 | Qui-quadrado | 2,0 × 10⁻⁵ | 0,257 |
| Raça × grau (Unknown excluído) | 355 | Fisher | 1,0 × 10⁻⁴ | 0,251 |
| Idade × grau | 363 | Wilcoxon | 4,9 × 10⁻³ | — |
| log2(AFP+1) × sexo | 276 | Wilcoxon | 5,4 × 10⁻³ | — |
| Fibrose Ishak × estágio | 199 | Wilcoxon | 6,7 × 10⁻³ | — |
| Vírus × grau (Unknown excluído) | 161 | Qui-quadrado | 1,1 × 10⁻² | 0,236 |
| AFP ≥400 × estágio | 263 | Qui-quadrado | 1,4 × 10⁻² | 0,152 |
| Fibrose avançada × estágio | 199 | Qui-quadrado | 1,4 × 10⁻² | 0,174 |
| Vírus × estágio (Unknown excluído) | 156 | Qui-quadrado | 1,4 × 10⁻² | 0,233 |

### Associações não significativas

- **Sexo × estágio** (p = 0,163) e **sexo × grau** (p = 0,404).
- **Idade × estágio** (p = 0,035; BH = 0,060) e **idade × sexo** (p = 0,089).
- **Invasão vascular × estágio** (p = 0,079) e **× grau** (p = 0,363).
- **Fibrose × grau** (p = 0,504).
- **Histologia HCC vs. outro × estágio** (p = 1,0).

**Interpretação:**

- **Child-Pugh** e **estágio** estão fortemente associados: pacientes com
  doença mais avançada têm pior função hepática.
- **AFP** se associa mais ao **grau** (ρ = 0,292) do que ao estágio
  (ρ = 0,162), consistente com produção de AFP em tumores menos
  diferenciados.
- **Fibrose** se associa ao estágio, mas não ao grau.
- **Sexo** não se associa a estágio/grau nesta coorte.
- A associação **raça × grau** deve ser interpretada com cautela: raça é uma
  construção social e pode refletir ancestralidade, acesso ao cuidado e
  diferenças de recrutamento.

---

## 8.4 AFP

![AFP por estágio e grau](resultados/figuras/15_AFP_por_estagio_grau.png)

*Figura 15. log2(AFP + 1) por estágio agrupado e por grau agrupado (violino + boxplot).*

Medianas de AFP (ng/mL):

| Subgrupo | n | Mediana | Q1–Q3 |
|---|---:|---:|---:|
| Estágio I–II | 211 | 12,0 | 4,0–188,5 |
| Estágio III–IV | 52 | 34,5 | 5,0–2.628,8 |
| Grau G1–G2 | 162 | 7,0 | 4,0–47,8 |
| Grau G3–G4 | 113 | 48,0 | 9,0–1.902,0 |
| Feminino | 91 | 28,0 | 6,0–805,0 |
| Masculino | 185 | 11,0 | 4,0–120,0 |

A **média** de AFP é muito inflada por outliers (por exemplo, média de
14.622 ng/mL no estágio I–II), por isso a interpretação usa **mediana** e
**log2(AFP+1)**.

---

## 8.5 Sobrevida global (Kaplan-Meier)

Foram usados 363 tumores com `OS_time > 0` e `OS` informado; 128 eventos.

| Variável | n | Eventos | log-rank p |
|---|---:|---:|---:|
| Estágio (I–II vs. III–IV) | 339 | 114 | **5,18 × 10⁻⁷** |
| Grau (G1–G2 vs. G3–G4) | 358 | 124 | 0,682 |
| Sexo | 363 | 128 | 0,288 |
| Idade (≥65 vs. <65) | 363 | 128 | 0,244 |
| AFP (≥400 vs. <400) | 274 | 82 | 0,941 |
| Fibrose avançada (Ishak ≥5) | 207 | 62 | 0,428 |

![KM estágio](resultados/figuras/17_KM_OS_estagio.png)

*Figura 17. Sobrevida global por estágio.*

![KM grau](resultados/figuras/18_KM_OS_grau.png)

*Figura 18. Sobrevida global por grau.*

![KM sexo](resultados/figuras/19_KM_OS_sexo.png)

*Figura 19. Sobrevida global por sexo.*

![KM idade](resultados/figuras/20_KM_OS_idade.png)

*Figura 20. Sobrevida global por idade.*

![KM AFP](resultados/figuras/21_KM_OS_AFP.png)

*Figura 21. Sobrevida global por AFP.*

![KM fibrose](resultados/figuras/22_KM_OS_fibrose.png)

*Figura 22. Sobrevida global por fibrose.*

![Painel KM](resultados/figuras/23_KM_OS_painel.png)

*Figura 23. Painel comparativo das curvas de estágio, grau, sexo e idade.*

**Interpretação:** apenas o **estágio** separou claramente as curvas de
sobrevida. A ausência de significância para grau, AFP e fibrose pode refletir
tanto ausência de efeito prognóstico independente quanto perda de poder por
missingness e heterogeneidade.

---

## 8.6 Cox univariado

A tabela completa está em `resultados/tabelas/17_cox_univariado.csv`.

| Variável | HR (IC 95%) | p |
|---|---|---:|
| Estágio III–IV vs. I–II | 2,52 (1,73–3,66) | 1,2 × 10⁻⁶ |
| Grau G3–G4 vs. G1–G2 | 1,08 (0,75–1,55) | 0,682 |
| Sexo masculino vs. feminino | 0,82 (0,58–1,18) | 0,290 |
| Idade (por ano) | 1,012 (0,998–1,026) | 0,094 |
| AFP ≥400 vs. <400 | 1,02 (0,62–1,68) | 0,939 |
| Fibrose avançada vs. não | 0,80 (0,46–1,39) | 0,429 |
| Vírus HBV vs. nenhum/outro | 0,44 (0,25–0,76) | 3,1 × 10⁻³ |
| Vírus HCV vs. nenhum/outro | 0,57 (0,24–1,34) | 0,196 |
| Child-Pugh B–C vs. A | 2,36 (1,67–3,35) | 1,5 × 10⁻⁶ |

**Atenção:** o HR de HBV <1 no modelo univariado pode refletir confundimento
por estágio, função hepática e vigilância (diagnóstico mais precoce em
pacientes com HBV), e não um efeito protetor causal. O termo HBV+HCV não foi
estimável (ausência de eventos suficientes em uma das categorias). No modelo
multivariado, que não inclui vírus, o estágio permanece dominante.

---

## 8.7 Cox multivariado

![Forest plot Cox](resultados/figuras/24_forest_cox_multivariado.png)

*Figura 24. Hazard ratios do Cox multivariado (Modelos A e B).*

### Modelo A: estágio + grau + sexo + idade

n = 337; 112 eventos.

| Termo | HR (IC 95%) | p |
|---|---|---:|
| Estágio III–IV vs. I–II | 2,61 (1,79–3,81) | 6,0 × 10⁻⁷ |
| Grau G3–G4 vs. G1–G2 | 1,16 (0,79–1,71) | 0,450 |
| Sexo masculino vs. feminino | 0,88 (0,59–1,29) | 0,508 |
| Idade (por ano) | 1,011 (0,996–1,026) | 0,155 |

### Modelo B: Modelo A + AFP + fibrose

n = 179; 52 eventos.

| Termo | HR (IC 95%) | p |
|---|---|---:|
| Estágio III–IV vs. I–II | 2,23 (1,24–4,03) | 7,7 × 10⁻³ |
| Grau G3–G4 vs. G1–G2 | 1,47 (0,84–2,58) | 0,176 |
| Sexo masculino vs. feminino | 0,99 (0,54–1,84) | 0,983 |
| Idade (por ano) | 1,031 (1,006–1,058) | 0,0165 |
| AFP ≥400 vs. <400 | 1,04 (0,51–2,12) | 0,907 |
| Fibrose avançada vs. não | 1,16 (0,62–2,17) | 0,645 |

**Interpretação:** o estágio é o preditor independente mais forte. O grau perde
significância após ajuste pelo estágio. No Modelo B, a idade aparece como fator
independente (HR ≈ 1,03 por ano), mas o modelo tem apenas 52 eventos e n = 179;
o achado é frágil. AFP e fibrose não foram independentes após ajuste.

---

## 8.8 Riscos proporcionais (`cox.zph`)

| Modelo | Termo | χ² | p |
|---|---|---:|---:|
| A | Estágio | 0,53 | 0,465 |
| A | Grau | 0,43 | 0,514 |
| A | Sexo | 2,70 | 0,100 |
| A | Idade | 0,39 | 0,532 |
| A | **GLOBAL** | 4,15 | **0,386** |
| B | Estágio | 2,17 | 0,141 |
| B | Grau | 0,86 | 0,355 |
| B | Sexo | 1,39 | 0,238 |
| B | Idade | 0,09 | 0,758 |
| B | AFP | 3,88 | 0,049 |
| B | Fibrose | 2,07 | 0,150 |
| B | **GLOBAL** | 9,38 | **0,153** |

Não há evidência global de violação da suposição de riscos proporcionais. AFP
no Modelo B está no limite (p = 0,049), o que merece cautela.

---

## 8.9 Correlações clínicas

![Correlações clínicas](resultados/figuras/25_correlacao_clinica.png)

*Figura 25. Correlações de Spearman entre idade, AFP, fibrose, estágio e grau.*

| Par | ρ | p | n |
|---|---:|---:|---:|
| AFP × grau | 0,292 | 8,6 × 10⁻⁷ | 275 |
| AFP × estágio | 0,162 | 0,0083 | 263 |
| Fibrose × estágio | −0,205 | 0,0038 | 199 |
| Idade × grau | −0,131 | 0,0126 | 363 |
| Idade × estágio | −0,122 | 0,0235 | 344 |
| Idade × fibrose | −0,176 | 0,0109 | 209 |
| Fibrose × grau | 0,083 | 0,231 | 208 |
| Estágio × grau | −0,009 | 0,875 | 343 |

**Interpretação:** AFP e grau compartilham parte da variação; fibrose e estágio
são inversamente correlacionados (possivelmente porque pacientes com cirrose
avançada são diagnosticados em estágios diferentes ou têm maior mortalidade
competitiva). A correlação estágio × grau é praticamente nula, sugerindo que
eles capturam dimensões parcialmente distintas.

---

## 8.10 Limitações

1. **Missingness** elevado em invasão vascular, sorologia, fibrose e AFP.
2. **Viés de seleção:** a coorte TCGA é cirúrgica, majoritariamente estágio
   I–II; os resultados não se generalizam para todos os HCC.
3. **Dados observacionais:** associações não implicam causalidade.
4. **Múltiplos testes:** usamos BH, mas as análises são exploratórias.
5. **Categorizações:** estágio I–II vs. III–IV e grau G1–G2 vs. G3–G4 são
   escolhas pragmáticas; outros cortes podem mudar o poder.
6. **AFP:** distribuição muito assimétrica; mediana/log2 são mais adequados.
7. **Sobrevida:** não ajustamos para tratamento, performance status, função
   hepática detalhada ou causas competitivas de morte.
8. **Raça/etnia:** categorias pequenas e heterogêneas; interpretar com cautela.

---

## 8.11 Resumo

- O **estágio** é o fator prognóstico dominante na coorte TCGA-LIHC.
- **AFP** se associa a **grau**, não a estágio de forma independente.
- **Child-Pugh** e **fibrose** se associam ao estágio, refletindo função
  hepática e doença de base.
- **Sexo** não se associou a estágio, grau ou sobrevida nesta coorte.
- A **missingness** limita as conclusões sobre AFP, fibrose, vírus e invasão
  vascular.
- Os modelos de Cox não mostram violação global de riscos proporcionais.

---

# 9. Auditorias, benchmarks e compliance

> Arquivo: `docs/09_auditorias_benchmarks.md`.
> Scripts: `scripts/10_auditoria_dados_e_volcano.R`, `scripts/11_benchmarks.R`,
> `scripts/12_compliance.R`.
> Relatórios: `auditorias/`.

---

## 9.1 Auditoria de dados

Relatório completo: [`auditorias/AUDITORIA_DADOS.md`](docs/auditorias/AUDITORIA_DADOS.md).

Foram verificados:

- número e duplicidade de amostras;
- IDs TCGA/GTEx e exclusão das recidivas `-02`;
- dimensões, NAs, infinitos e duplicatas da matriz;
- alinhamento entre matriz e metadados;
- sexo, idade, estágio, grau e sobrevida;
- consistência entre variáveis derivadas;
- confundimento estudo × condição;
- existência dos arquivos essenciais.

**Resultado:** 22 verificações; todas OK, exceto o **alerta** de confundimento
estrutural (Cramér's V = 1). Não há NAs, infinitos, genes duplicados ou
amostras dessincronizadas na matriz.

**Missingness clínico (TCGA tumor):** invasão vascular 70,7%, sorologia viral
55,8%, fibrose 43,1%, Child-Pugh 35,8%, AFP 25,2%, estágio 6,5%, grau 1,4%,
idade 0,3%.

---

## 9.2 Auditoria do `Volcano_12`

Relatório completo:
[`auditorias/AUDITORIA_VOLCANO_12.md`](docs/auditorias/AUDITORIA_VOLCANO_12.md).

**Queixa:** a figura parecia vazia.

**Diagnóstico:** a imagem original não estava tecnicamente em branco
(≈3,65% de pixels não brancos), mas o ComBat sem covariável deixou apenas
**1 DEG** (LINC01419). Os outros 12.400 genes ficaram amontoados perto da
origem, tornando o painel visualmente vazio.

**Causa raiz:** lote (TCGA/GTEx) e condição (tumor/normal) são perfeitamente
confundidos; o ComBat removeu também a diferença biológica. A tentativa com
`mod = ~ grupo` é recusada pelo próprio ComBat.

**Correção:** a figura passou a ter **dois painéis** — antes (3.388 DEGs) e
depois (1 DEG) —, com escala livre, rótulo de LINC01419 e subtítulo com as
contagens. A versão corrigida tem ≈9,48% de pixels não brancos.

---

## 9.3 Benchmarks computacionais

Relatório completo: [`auditorias/BENCHMARKS.md`](docs/auditorias/BENCHMARKS.md).

| Etapa | Tempo | Observação |
|---|---:|---|
| Carga de metadados + matriz | 1,5 s | 12.401 genes × 529 amostras |
| Memória da matriz | 53,3 MB | matriz numérica |
| UMAP (2.000 genes) | 4,34 s | `n_neighbors=30`, `min_dist=0,25` |
| limma com covariáveis | 0,69 s | 12.401 genes × 478 amostras |
| limma sem covariáveis | 0,26 s | mesmo contraste |
| limma com `trend = FALSE` | 0,47 s | sensibilidade |

**Conclusão:** o custo dominante é a leitura da matriz e a geração de figuras;
o `limma` é rápido para esta escala.

---

## 9.4 Benchmarks de concordância estatística

| Comparação | Métrica | Valor |
|---|---|---:|
| Com vs. sem covariáveis | DEGs com / sem | 3.388 / 3.316 |
| Com vs. sem covariáveis | Jaccard | 0,948 |
| Com vs. sem covariáveis | ρ de logFC | 0,999 |
| `trend = TRUE` vs. `FALSE` | DEGs | 3.388 / 3.389 |
| `trend = TRUE` vs. `FALSE` | Jaccard | 1,000 |
| FDR apenas vs. FDR + \|logFC\| | genes | 10.109 / 3.388 |
| FDR apenas vs. FDR + \|logFC\| | Jaccard | 0,335 |
| Pareado vs. não pareado | DEGs | 3.710 / 4.130 |
| Pareado vs. não pareado | Jaccard | 0,764 |
| Pareado vs. não pareado | ρ de logFC | 0,969 |
| ComBat antes vs. depois | DEGs | 3.388 / 1 |
| ComBat antes vs. depois | Jaccard | 0,000 |
| ComBat antes vs. depois | ρ de logFC | 0,647 |

**Interpretação:**

- O ajuste por sexo e idade muda pouco o ranking global (ρ = 0,999), mas altera
  ~2% dos DEGs (Jaccard = 0,948). Reportar o modelo ajustado é mais defensável.
- `trend = TRUE` e `FALSE` produzem essencialmente o mesmo conjunto.
- O corte de \|logFC\| > 1 reduz de 10.109 para 3.388 genes; é uma escolha
  conservadora, não obrigatória.
- A análise pareada e a não pareada concordam substancialmente (Jaccard = 0,764;
  ρ = 0,969), reforçando o sinal tumor vs. adjacente.
- O ComBat sem covariável zera a concordância de DEGs (Jaccard = 0), confirmando
  a sobrecorreção.

---

## 9.5 Compliance

Relatório completo:
[`auditorias/COMPLIANCE.md`](docs/auditorias/COMPLIANCE.md).

Critérios:

- arquivos essenciais presentes;
- **nenhum PDF** em `resultados/figuras`;
- PNGs válidos e não vazios;
- `Volcano_12` com conteúdo visual mínimo;
- CSVs legíveis e colunas críticas sem NA;
- scripts R parseáveis;
- sementes de reprodutibilidade;
- imagens dos docs e do HTML existentes;
- docs sem menções a PDF ou a métodos lineares de redução de dimensionalidade;
- metadados com 529 amostras e sem recidivas.

O resultado detalhado é gerado automaticamente por
`scripts/12_compliance.R` e salvo em `auditorias/COMPLIANCE.md`.

---

## 9.6 Reprodutibilidade

- Sementes: `set.seed(2026)` em UMAP, volcano e fenótipos.
- `uwot::umap` com `n_threads = 1` e `n_sgd_threads = 1` para maior
  determinismo.
- Versões de R e pacotes em `docs/sessionInfo.txt`.
- Entrada/saída documentadas em `docs/07_reproducibilidade.md`.
- `run_all.R` executa 00 → 01 → 02 → 03 → 04 → 05 → 10 → 11 → 12.

---

# 10. Outros grupos dentro dos dados

> Arquivo: `docs/10_grupos_alternativos.md`.
> Script: `scripts/06_grupos_alternativos.R`.
> Saídas: `resultados/tabelas/21–22_*.csv`,
> `resultados/figuras/26–27_*.png`.

---

## 10.1 Resposta curta

Sim. Além de **tumor primário**, **tecido adjacente** e **fígado normal GTEx**,
os dados contêm:

1. **Tumores recidivados** (`Recurrent Tumor`, n = 2);
2. **subtipos histológicos** dentro do TCGA:
   - carcinoma hepatocelular clássico (n = 359);
   - carcinoma hepatocolangiocarcinoma misto (n = 7);
   - carcinoma fibrolamelar (n = 3);
3. **grupos de etiologia viral** (HBV, HCV, HBV+HCV, nenhum/outro, unknown);
4. **centros TCGA (TSS)**, que funcionam como grupos técnicos de lote;
5. dois grupos de tecido normal: **TCGA adjacente** e **GTEx normal**.

O grupo mais “escondido” é o de **recidivas**: ele aparece no arquivo bruto,
mas foi excluído da análise principal. Os subtipos histológicos e os grupos
virais são grupos biológicos reais, mas pequenos ou com missingness elevado.

---

## 10.2 Inventário

![Grupos nos dados](resultados/figuras/26_grupos_dados.png)

*Figura 26. (A) Tipos de amostra, incluindo as 2 recidivas. (B) Subtipos
histológicos no TCGA tumor (escala log10).*

| Categoria | Grupo | n | % |
|---|---|---:|---:|
| Tipo de amostra | Tumor primário | 369 | 69,5 |
| Tipo de amostra | Tecido adjacente | 50 | 9,4 |
| Tipo de amostra | GTEx normal | 110 | 20,7 |
| Tipo de amostra | **Recidiva** | **2** | **0,4** |
| Histologia (TCGA tumor) | HCC clássico | 359 | 97,3 |
| Histologia (TCGA tumor) | Hepatocolangiocarcinoma misto | 7 | 1,9 |
| Histologia (TCGA tumor) | Fibrolamelar | 3 | 0,8 |
| Sorologia viral (TCGA tumor) | HBV | 59 | 16,0 |
| Sorologia viral (TCGA tumor) | HCV | 19 | 5,1 |
| Sorologia viral (TCGA tumor) | HBV+HCV | 85 | 23,0 |
| Sorologia viral (TCGA tumor) | Unknown | 206 | 55,8 |

Além disso, há **38 centros TCGA (TSS)**; o maior é `DD` (n = 144; 39,0% dos
tumores), seguido de `CC` e `G3` (n = 32 cada). Isso é relevante como possível
fonte de efeito de lote técnico, embora não tenha sido usado como covariável
nas análises principais.

---

## 10.3 Recidivas

As 2 amostras `TCGA-DD-AACA-02` e `TCGA-ZS-A9CF-02` são **tumores
recidivados**, não primários. Foram excluídas porque:

- o desenho principal é tumor primário vs. tecido adjacente/normal;
- n = 2 não permite análise estatística robusta;
- recidivas podem ter biologia distinta (resistência, evolução clonal).

Elas estão mantidas em `metadata_amostras.tsv` para auditoria, mas não em
`metadata_analise.tsv`.

---

## 10.4 Subtipos histológicos

O TCGA-LIHC contém, além do HCC clássico:

- **Hepatocolangiocarcinoma misto** (n = 7): tumor com componente
  hepatocelular e colangiocelular;
- **Carcinoma fibrolamelar** (n = 3): variante rara, geralmente em pacientes
  jovens, com biologia distinta.

Com n = 7 e n = 3, **não é adequado rodar expressão diferencial inferencial**
contra o HCC clássico. A análise é descritiva.

![Marcadores por histologia](resultados/figuras/27_histologia_marcadores.png)

*Figura 27. Expressão de marcadores por subtipo histológico. Cada ponto é uma
amostra; a barra preta é a mediana. Interpretar com cautela devido ao n
pequeno.*

Medianas de expressão (log2 TPM + 0,001):

| Gene | HCC clássico (n = 359) | Misto (n = 7) | Fibrolamelar (n = 3) |
|---|---:|---:|---:|
| GPC3 | 7,28 | 3,10 | 4,58 |
| AFP | 3,53 | 2,93 | 1,95 |
| KRT19 | −0,28 | **7,24** | −0,81 |
| EPCAM | −0,23 | **6,53** | 2,10 |
| BIRC5 | 2,77 | 3,79 | 0,62 |
| VIM | 6,98 | 8,48 | 8,94 |

**Interpretação descritiva:**

- O grupo **misto** tem expressão muito alta de **KRT19** e **EPCAM**,
  consistente com o componente colangiocelular/ductal.
- O grupo **fibrolamelar** tem AFP e GPC3 mais baixos e VIM alto, mas n = 3
  impede conclusões.
- O HCC clássico tem a maior mediana de GPC3, coerente com o marcador
  hepatocelular.
- Esses achados são **hipóteses exploratórias**, não resultados confirmatórios.

---

## 10.5 Grupos virais

A sorologia viral disponível no TCGA mostra:

- **HBV+HCV:** 85 tumores (23,0%);
- **HBV:** 59 (16,0%);
- **HCV:** 19 (5,1%);
- **Unknown:** 206 (55,8%).

O grupo `Unknown` é grande e não deve ser interpretado como “sem infecção”.
Nas análises de fenótipos, a categoria `Unknown` foi excluída dos testes de
associação e do Cox univariado. A etiologia viral pode ser usada como
estratificação em análises futuras, mas o missingness limita o poder.

---

## 10.6 Centros TCGA (TSS)

Os 38 centros TCGA são um grupo **técnico**, não biológico. O centro `DD`
contribui com 39% dos tumores; centros com n = 1–3 são muitos. Isso pode gerar
efeito de lote. Não aplicamos ComBat por TSS porque:

- o efeito de lote principal neste estudo é TCGA vs. GTEx;
- o número de amostras por centro é muito desigual;
- o desenho não foi pré-especificado para correção por TSS.

Ainda assim, o TSS é um candidato natural para análises de sensibilidade
futuras (por exemplo, incluí-lo como efeito aleatório ou corrigir apenas
centros com n ≥ 10).

---

## 10.7 O que isso muda na interpretação?

- A análise principal continua válida para **tumor primário vs. adjacente vs.
  GTEx normal**.
- As recidivas não entram nas comparações principais.
- Os subtipos histológicos e virais são **grupos reais**, mas pequenos ou com
  missingness; não foram usados para inferência principal.
- A existência de 38 TSS reforça a necessidade de cautela com efeitos de lote
  e de reportar a origem das amostras.
- Uma análise futura interessante é **estratificar por etiologia viral** ou
  por **subtipo histológico** em coortes maiores.

---

## 10.8 Resumo

| Grupo | n | Pode ser analisado? |
|---|---:|---|
| Tumor primário | 369 | Sim (principal) |
| Tecido adjacente | 50 | Sim (principal) |
| GTEx normal | 110 | Sim (principal) |
| Recidiva | 2 | Não (descritivo/apenas auditoria) |
| HCC clássico | 359 | Sim |
| Misto | 7 | Apenas descritivo |
| Fibrolamelar | 3 | Apenas descritivo |
| HBV / HCV / HBV+HCV | 59 / 19 / 85 | Possível, com cautela |
| Unknown viral | 206 | Excluído dos testes |
| TSS | 38 centros | Efeito técnico; não usado |

As análises separadas de cada um desses grupos estão em
[`docs/11_grupos_separados.md`](docs/11_grupos_separados.md): clínica das
recidivas, expressão no painel LA, subtipos histológicos, grupos virais e
sinal de lote por TSS no UMAP.

---

# 11. Análises separadas dos grupos alternativos

> Arquivo: `docs/11_grupos_separados.md`.
> Script: `scripts/07_grupos_separados.R`.
> Saídas: `resultados/tabelas/23–35_*.csv`,
> `resultados/figuras/28–35_*.png` e
> `resultados/qc/resumo_grupos_separados.txt`.

---

## 11.1 Princípio geral

A análise principal continua sendo **tumor primário × tecido adjacente ×
fígado normal GTEx**. Os grupos analisados neste capítulo são mantidos
**separados**, sem serem fundidos aos três grupos principais.

Como as recidivas têm n = 2, os subtipos histológicos raros têm n = 7 e n = 3,
e a sorologia viral tem 55,8% de `Unknown`, essas análises são
**descritivas/exploratórias**. Elas não substituem testes confirmatórios e não
entram no desenho principal.

| Grupo | n | Tratamento neste capítulo |
|---|---:|---|
| Tumores recidivados | 2 | Descritivo; expressão no painel LA |
| HCC clássico | 359 | Referência dos subtipos |
| Hepatocolangiocarcinoma misto | 7 | Descritivo/exploratório |
| Carcinoma fibrolamelar | 3 | Descritivo/exploratório |
| HBV | 59 | Exploratório, com `Unknown` separado |
| HCV | 19 | Exploratório, com `Unknown` separado |
| HBV+HCV | 85 | Exploratório, com `Unknown` separado |
| Unknown viral | 206 | Mostrado, mas excluído dos testes virais |
| Centros TCGA (TSS) | 38 | Grupo técnico de lote |

---

## 11.2 Tumores recidivados (n = 2)

### 11.2.1 Clínica

As duas amostras `Recurrent Tumor` são:

| Amostra | Idade | Sexo | Estágio | Grau | AFP | Fibrose Ishak | Vírus | Óbito | TSS |
|---|---:|---|---|---|---:|---:|---|---|---|
| `TCGA-DD-AACA-02` | 65 | masculino | I | G3 | 7 | 6 | Unknown | não | DD |
| `TCGA-ZS-A9CF-02` | 64 | masculino | II | G2 | NA | 0 | HBV+HCV | não | ZS |

Comparação descritiva com os 369 tumores primários:

| Variável | Recidiva 1 | Recidiva 2 | Primários (n = 369) |
|---|---|---|---|
| Idade (anos) | 65,0 | 64,0 | 61,0 (51,8–69,0) |
| Sexo masculino | sim | sim | 249 (67,5%) |
| Estágio I–II | sim | sim | 255 (73,9%) |
| Grau G3–G4 | sim | não | 132 (36,3%) |
| AFP (ng/mL) | 7,0 | NA | 15,5 (4,0–270,5) |
| Fibrose Ishak | 6,0 | 0,0 | 3,0 |
| Vírus detectável | não (Unknown) | sim (HBV+HCV) | 163 (44,2%) |
| Óbito | não | não | 128 (34,7%) |

Tabelas completas: `resultados/tabelas/23_recidivas_clinico.csv` e
`resultados/tabelas/24_recidivas_comparacao.csv`.

As duas recidivas são homens, com estágios iniciais (I e II), graus G3/G2 e
ambos vivos no último seguimento (2.301 e 2.412 dias). A idade e o estágio não
são atípicos em relação aos primários, mas **n = 2 não permite qualquer
inferência**.

### 11.2.2 Expressão no painel LA

A matriz completa (`expressao_log2tpm_filtrada.tsv.gz`) contém 529 amostras e
**não inclui as recidivas `-02`**. Por isso, a expressão das recidivas foi
analisada no **painel LA** (`expressao_LA_pathway.tsv.gz`, 213 genes), que
contém as 531 amostras.

Para cada gene, calcularam-se o z-score e o percentil da recidiva em relação
aos 369 tumores primários. A figura mostra os 10 genes com maior |z| em pelo
menos uma recidiva.

![Expressão das recidivas no painel LA](resultados/figuras/28_recidivas_expressao.png)

*Figura 28. Top 10 genes do painel LA com maior |z| nas recidivas. Violino =
distribuição dos 369 tumores primários; pontos = recidivas.*

Principais extremos (tabela completa em
`resultados/tabelas/25_recidivas_expressao.csv`):

| Recidiva | Gene | z | Percentil entre primários |
|---|---|---:|---:|
| `TCGA-DD-AACA-02` | LBP | −2,85 | 1,9 |
| `TCGA-DD-AACA-02` | PRKCA | 2,73 | 100,0 |
| `TCGA-DD-AACA-02` | NFATC3 | 2,52 | 99,2 |
| `TCGA-DD-AACA-02` | MAP2K6 | −2,41 | 0,8 |
| `TCGA-DD-AACA-02` | NCF4 | −2,34 | 0,8 |
| `TCGA-ZS-A9CF-02` | CAMK2D | 3,07 | 100,0 |
| `TCGA-ZS-A9CF-02` | PIK3CA | 2,88 | 99,7 |
| `TCGA-ZS-A9CF-02` | NFATC3 | 2,64 | 99,7 |
| `TCGA-ZS-A9CF-02` | ATF6 | 2,63 | 99,2 |
| `TCGA-ZS-A9CF-02` | PPP3CA | 2,38 | 99,5 |

Os genes extremos são do painel LA (inflamação, lipídios e sinalização), **não
são marcadores canônicos de HCC**. Alguns extremos são compartilhados
(`NFATC3` alto nas duas recidivas) e outros são específicos de uma amostra
(`LBP` muito baixo na primeira; `CAMK2D`/`PIK3CA` altos na segunda). Com n = 2,
isso é apenas uma descrição.

### 11.2.3 Similaridade global com os grupos

A correlação de Spearman foi calculada entre cada recidiva e todas as amostras
de cada grupo, usando os 213 genes do painel LA:

| Recidiva | Grupo | n | ρ mediana | ρ máximo |
|---|---|---:|---:|---:|
| `TCGA-DD-AACA-02` | Primário | 369 | 0,885 | 0,966 |
| `TCGA-DD-AACA-02` | Adjacente | 50 | 0,883 | 0,912 |
| `TCGA-DD-AACA-02` | Normal | 110 | 0,871 | 0,902 |
| `TCGA-ZS-A9CF-02` | Primário | 369 | 0,877 | 0,975 |
| `TCGA-ZS-A9CF-02` | Adjacente | 50 | 0,925 | 0,945 |
| `TCGA-ZS-A9CF-02` | Normal | 110 | 0,903 | 0,932 |

![Correlação global das recidivas](resultados/figuras/29_recidivas_correlacao.png)

*Figura 29. Distribuição das correlações de cada recidiva com primários,
adjacentes e normais.*

As correlações medianas são altas em todos os grupos (> 0,87), como esperado
para tecido hepático. A segunda recidiva é mais próxima do tecido adjacente
(ρ = 0,925) do que dos primários (ρ = 0,877), mas isso não deve ser
interpretado como evidência biológica: são duas amostras e a diferença pode
refletir composição tecidual ou ruído.

Tabela completa: `resultados/tabelas/26_recidivas_correlacao.csv`.

**Limitações das recidivas:** n = 2; apenas o painel LA; ausência de metadados
de tratamento prévio; as recidivas não foram emparelhadas com o tumor primário
do mesmo paciente na matriz disponível.

---

## 11.3 Subtipos histológicos

### 11.3.1 Composição clínica

| Subtipo | n | Idade mediana (Q1–Q3) | Masculino | Estágio III–IV | Grau G3–G4 | AFP mediana (Q1–Q3) | Óbitos |
|---|---:|---:|---:|---:|---:|---:|---:|
| HCC clássico | 359 | 61 (52–69) | 67,7% | 26,3% | 36,7% | 16 (4–267) | 127 |
| Misto | 7 | 62 (61–70) | 57,1% | 14,3% | 14,3% | 14 (3–290) | 1 |
| Fibrolamelar | 3 | 25 (22,5–29,5) | 66,7% | 33,3% | 33,3% | 2,5 (2,25–2,75) | 0 |

Tabela completa: `resultados/tabelas/27_histologia_clinico.csv`.

O achado clínico mais evidente é a **idade muito menor no carcinoma
fibrolamelar** (mediana 25 anos vs. 61–62 anos nos outros grupos) e a **AFP
baixa** nesse subtipo. O grupo misto tem apenas 7 casos e o fibrolamelar 3;
as proporções de estágio, grau e óbito são instáveis.

![Idade e AFP por subtipo histológico](resultados/figuras/30_histologia_clinico.png)

*Figura 30. Distribuição de idade e de log2(AFP + 1) por subtipo histológico.*

### 11.3.2 Marcadores

A figura 27 (em `docs/10_grupos_alternativos.md`) mostra a expressão de GPC3,
AFP, KRT19, EPCAM, BIRC5 e VIM por subtipo. A tabela de medianas está em
`resultados/tabelas/22_histologia_marcadores.csv`.

Testes de Kruskal-Wallis entre os três subtipos, com BH para os nove marcadores
(`resultados/tabelas/28_histologia_testes.csv`):

| Gene | p Kruskal | p BH | Direção descritiva |
|---|---:|---:|---|
| KRT19 | 0,00061 | **0,0018** | muito alto no misto |
| EPCAM | 0,0047 | **0,0106** | muito alto no misto |
| VIM | 0,000082 | **0,00074** | alto no misto e no fibrolamelar |
| COL1A1 | 0,00037 | **0,0017** | alto no misto e no fibrolamelar |
| BIRC5 | 0,0398 | 0,0598 | sem diferença após BH |
| SPP1 | 0,0312 | 0,0561 | sem diferença após BH |
| MMP9 | 0,0530 | 0,0681 | sem diferença após BH |
| GPC3 | 0,334 | 0,334 | sem diferença |
| AFP | 0,135 | 0,152 | sem diferença |

Os destaques biológicos são coerentes com a histologia:

- o **misto** tem KRT19 e EPCAM muito elevados, consistente com componente
  ductal/colangiocelular;
- **VIM** e **COL1A1** são mais altos nos subtipos não clássicos, sugerindo
  maior remodelamento/estroma;
- GPC3 e AFP não diferenciam os grupos de forma robusta.

**Limitações:** n = 7 (misto) e n = 3 (fibrolamelar). Os p-valores são
exploratórios; não há poder para expressão diferencial genômica confiável.
Não foi feita inferência transcriptoma-wide por subtipo.

---

## 11.4 Grupos virais

### 11.4.1 Composição clínica

| Grupo viral | n | Idade mediana | Masculino | Estágio III–IV | Grau G3–G4 | AFP mediana | Fibrose mediana | Óbitos |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HBV | 59 | 64 | 61,0% | 30,4% | 42,4% | 11 | 0 | 22 |
| HCV | 19 | 64 | 84,2% | 11,8% | 11,1% | 11 | 1 | 6 |
| HBV+HCV | 85 | 59 | 71,8% | 45,8% | 23,8% | 7 | 1 | 42 |
| Unknown | 206 | 61 | 66,0% | 17,5% | 41,9% | 24,5 | 3 | 58 |

Tabela completa: `resultados/tabelas/29_viral_clinico.csv`.

### 11.4.2 Testes e missingness

Testes entre os três grupos conhecidos (`HBV`, `HCV`, `HBV+HCV`) e entre
**conhecido vs. `Unknown`** estão em `resultados/tabelas/30_viral_testes.csv`.
Após BH:

| Comparação | Variável | p BH | Interpretação |
|---|---|---:|---|
| HBV vs HCV vs HBV+HCV | Fibrose Ishak | **0,0389** | diferença entre grupos conhecidos |
| HBV vs HCV vs HBV+HCV | Estágio | **0,0346** | diferença entre grupos conhecidos |
| HBV vs HCV vs HBV+HCV | Grau | **0,0346** | diferença entre grupos conhecidos |
| HBV vs HCV vs HBV+HCV | Idade | 0,376 | sem diferença |
| HBV vs HCV vs HBV+HCV | AFP | 0,976 | sem diferença |
| HBV vs HCV vs HBV+HCV | Sexo | 0,190 | sem diferença |
| Conhecido vs Unknown | Estágio | **0,00071** | `Unknown` tem menos III–IV |
| Conhecido vs Unknown | Grau | **0,0346** | `Unknown` tem mais G3–G4 |
| Conhecido vs Unknown | Fibrose | **0,0023** | `Unknown` tem mediana maior (3 vs. 0–1) |
| Conhecido vs Unknown | Idade | 0,784 | sem diferença |
| Conhecido vs Unknown | AFP | 0,190 | sem diferença |
| Conhecido vs Unknown | Sexo | 0,601 | sem diferença |

O ponto crítico é que **`Unknown` não é missing aleatório**: ele difere dos
grupos conhecidos em estágio, grau e fibrose. Portanto, análises virais
restritas aos casos conhecidos podem ser enviesadas. `Unknown` **não deve ser
interpretado como “sem infecção”**.

### 11.4.3 Marcadores por grupo viral

As medianas de GPC3, AFP, KRT19, EPCAM, BIRC5, VIM, SPP1, MMP9 e COL1A1 por
grupo viral estão em `resultados/tabelas/31_viral_marcadores.csv`.

Nenhum dos nove marcadores diferiu entre HBV, HCV e HBV+HCV após BH (menor
p BH = 0,566, para BIRC5). Ou seja, neste conjunto, **a etiologia viral não
se traduz em diferença clara de expressão desses marcadores**.

![AFP por grupo viral](resultados/figuras/31_viral_clinico.png)

*Figura 31. log2(AFP + 1) por grupo viral. `Unknown` é mostrado, mas não entra
nos testes.*

![Marcadores por grupo viral](resultados/figuras/32_viral_marcadores.png)

*Figura 32. Marcadores por grupo viral. Violino + boxplot; testes apenas entre
HBV, HCV e HBV+HCV.*

**Limitações:** 55,8% dos tumores têm sorologia `Unknown`; a definição de
“HBV+HCV” pode refletir coinfecção ou marcadores sorológicos; não há ajuste
para estágio, grau, tratamento ou TSS nesta análise descritiva.

---

## 11.5 Centros TCGA (TSS)

### 11.5.1 Distribuição

O TCGA-LIHC tem **38 centros (TSS)**. A distribuição é muito desigual:

| TSS | n | % dos tumores | Idade mediana | Estágio III–IV | Grau G3–G4 | AFP mediana |
|---|---:|---:|---:|---:|---:|---:|
| DD | 144 | 39,0% | 61 | 17,4% | 53,5% | 11 |
| CC | 32 | 8,7% | 54,5 | 81,3% | 21,9% | NA |
| G3 | 32 | 8,7% | 58,5 | 18,8% | 21,9% | 6,5 |
| BC | 21 | 5,7% | 69 | 50,0% | 45,0% | 248,5 |
| 2Y | 20 | 5,4% | 64 | 5,6% | 15,0% | 19 |
| ED | 15 | 4,1% | 53 | 46,7% | 53,3% | 145,5 |
| FV | 11 | 3,0% | 75 | 9,1% | 0,0% | 67 |
| Outros (31 TSS) | 94 | 25,5% | — | — | — | — |

Tabela completa com os 38 centros: `resultados/tabelas/32_tss_resumo.csv`.

Sete centros têm n ≥ 10 e foram usados nos testes: **DD, CC, G3, BC, 2Y, ED,
FV**. Os demais foram mantidos no resumo, mas não em testes específicos por
centro, devido ao n muito pequeno.

![Distribuição por TSS](resultados/figuras/33_tss_distribuicao.png)

*Figura 33. Número de tumores por TSS (top 14 + agregado “Outros”).*

### 11.5.2 Clínica e marcadores entre centros

Testes entre os sete centros com n ≥ 10
(`resultados/tabelas/33_tss_testes.csv`):

| Variável | Teste | n | p BH |
|---|---|---:|---:|
| Idade | Kruskal-Wallis | 275 | **4,2 × 10⁻⁵** |
| Sexo | Qui-quadrado/Fisher | 275 | **0,0476** |
| Estágio | Qui-quadrado/Fisher | 264 | **0,00030** |
| Grau | Qui-quadrado/Fisher | 272 | **0,00050** |
| AFP | Kruskal-Wallis | 213 | 0,0614 |
| Fibrose Ishak | Kruskal-Wallis | 154 | 0,0614 |
| KRT19 | Kruskal-Wallis | 275 | **2,8 × 10⁻⁵** |
| EPCAM | Kruskal-Wallis | 275 | **2,8 × 10⁻⁵** |
| BIRC5 | Kruskal-Wallis | 275 | **2,8 × 10⁻⁵** |
| COL1A1 | Kruskal-Wallis | 275 | **0,0018** |
| GPC3 | Kruskal-Wallis | 275 | 0,193 |
| VIM | Kruskal-Wallis | 275 | 0,116 |
| SPP1 | Kruskal-Wallis | 275 | 0,453 |
| MMP9 | Kruskal-Wallis | 275 | 0,099 |
| AFP (expressão) | Kruskal-Wallis | 275 | 0,0614 |

Os centros diferem em **idade, sexo, estágio, grau** e em marcadores como
**KRT19, EPCAM, BIRC5 e COL1A1**. Isso é compatível com efeito de centro/lote
e também com composição clínica desigual entre centros; as duas explicações
não são separáveis neste desenho.

![Medianas de marcadores por TSS](resultados/figuras/35_tss_marcadores_heatmap.png)

*Figura 35. Medianas de marcadores por TSS (z-score entre centros), para os
sete centros com n ≥ 10. A anotação superior mostra o n de cada centro.*

### 11.5.3 Sinal de lote no UMAP

As coordenadas UMAP já calculadas para as 529 amostras foram usadas para
verificar se o TSS se associa à posição no espaço de expressão. Foi comparado
o conjunto dos sete centros com n ≥ 10 contra “Outros” (ANOVA e
Kruskal-Wallis):

| Coordenada | R² do TSS | p ANOVA | p Kruskal-Wallis |
|---|---:|---:|---:|
| UMAP1 | 0,083 | 5,2 × 10⁻⁵ | 7,5 × 10⁻⁵ |
| UMAP2 | 0,072 | 3,3 × 10⁻⁴ | 2,6 × 10⁻⁴ |

O TSS explica cerca de **8,3% da variação de UMAP1** e **7,2% de UMAP2**
(p < 0,001). Não é a maior fonte de estrutura — grupos biológicos e
tumor/normal dominam —, mas é um sinal não desprezível.

![UMAP por TSS](resultados/figuras/34_tss_umap.png)

*Figura 34. UMAP dos tumores colorido por TSS. Centros com n ≥ 10 em cores;
demais agregados em “Outros”.*

**Limitações:** TSS não é randomizado; centros pequenos foram agregados; o
UMAP não é uma medida linear de expressão e o R² é descritivo. O TSS não foi
usado como covariável na análise principal; os resultados sugerem que uma
análise de sensibilidade com efeito de centro (ou correção restrita a centros
grandes) é justificável em trabalhos futuros.

---

## 11.6 Síntese

| Grupo | n | O que foi feito | Pode concluir? |
|---|---:|---|---|
| Recidivas | 2 | Clínica, expressão no painel LA e correlação global | Não; apenas descrição |
| HCC clássico | 359 | Referência clínica e de marcadores | Sim, como grupo principal |
| Misto | 7 | Clínica e marcadores | Apenas hipóteses |
| Fibrolamelar | 3 | Clínica e marcadores | Apenas hipóteses |
| HBV / HCV / HBV+HCV | 59 / 19 / 85 | Clínica, testes e marcadores | Exploratório; missingness informativa |
| Unknown viral | 206 | Mostrado separadamente | Não tratar como “sem infecção” |
| TSS | 38 | Distribuição, clínica, marcadores e UMAP | Efeito técnico plausível; não causal |

A conclusão geral é que **os grupos alternativos são reais e visíveis nos
dados**, mas nenhum deles tem tamanho ou qualidade de missingness suficiente
para substituir o desenho principal. O valor deste capítulo é documentar
heterogeneidade, levantar hipóteses e justificar cautela com efeitos de centro
e com a sorologia viral incompleta.

---

## 11.7 Reprodutibilidade

```bash
Rscript scripts/07_grupos_separados.R
```

O script lê `metadata_amostras.tsv`, `metadata_analise.tsv`,
`expressao_log2tpm_filtrada.tsv.gz`, `expressao_LA_pathway.tsv.gz` e
`resultados/tabelas/09_UMAP_coordenadas.csv`. As sementes são fixadas com
`set.seed(2026)`. O resumo textual é gravado em
`resultados/qc/resumo_grupos_separados.txt`.
