# 9. Auditorias, benchmarks e compliance

> Arquivo: `docs/09_auditorias_benchmarks.md`.
> Scripts: `scripts/10_auditoria_dados_e_volcano.R`, `scripts/11_benchmarks.R`,
> `scripts/12_compliance.R`.
> Relatórios: `auditorias/`.

---

## 9.1 Auditoria de dados

Relatório completo: [`auditorias/AUDITORIA_DADOS.md`](../auditorias/AUDITORIA_DADOS.md).

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
[`auditorias/AUDITORIA_VOLCANO_12.md`](../auditorias/AUDITORIA_VOLCANO_12.md).

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

Relatório completo: [`auditorias/BENCHMARKS.md`](../auditorias/BENCHMARKS.md).

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
[`auditorias/COMPLIANCE.md`](../auditorias/COMPLIANCE.md).

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
