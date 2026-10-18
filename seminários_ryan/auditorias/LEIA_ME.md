# Auditorias, benchmarks e compliance

| Arquivo | Conteúdo |
|---|---|
| `AUDITORIA_DADOS.md` | integridade dos dados, amostras, matriz, missingness e confundimento |
| `AUDITORIA_VOLCANO_12.md` | diagnóstico da figura `Volcano_12` e correção aplicada |
| `BENCHMARKS.md` | tempos de execução e concordância estatística |
| `COMPLIANCE.md` | checklist de conformidade dos artefatos |
| `auditoria_dados_resumo.csv` | tabela da auditoria de dados |
| `auditoria_volcano12_diagnostico.csv` | métricas da imagem e dos DEGs |
| `BENCHMARKS_tabela.csv` | tabela dos benchmarks |
| `COMPLIANCE_tabela.csv` | tabela do checklist |

## Como reproduzir

```bash
Rscript scripts/10_auditoria_dados_e_volcano.R
Rscript scripts/11_benchmarks.R
Rscript scripts/12_compliance.R
```

Os relatórios são regerados a partir dos dados e figuras atuais.
