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
