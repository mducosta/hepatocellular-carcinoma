# Auditorias, benchmarks e compliance

| Arquivo | Conteúdo |
|---|---|
| `auditoria_dados.md` | integridade dos dados, amostras, matriz, missingness e confundimento |
| `auditoria_volcano_ComBat.md` | diagnóstico da figura `volcano_12` e correção aplicada |
| `benchmarks_pipeline.md` | tempos de execução e concordância estatística |
| `compliance_pipeline.md` | checklist de conformidade dos artefatos |
| `auditoria_dados_resumo.csv` | tabela da auditoria de dados |
| `auditoria_volcano_ComBat_diagnostico.csv` | métricas da imagem e dos DEGs |
| `benchmarks_pipeline_tabela.csv` | tabela dos benchmarks |
| `compliance_pipeline_tabela.csv` | tabela do checklist |
| `auditoria_analises_avancadas.md` | auditoria crítica das covariáveis, ComBat e deconvolução |
| `auditoria_confundimento_associacoes.csv` | Cramér's V e p para pares de variáveis |
| `auditoria_fluxo_casos.csv` | fluxo de casos completos |
| `auditoria_atenuacao_sinais.csv` | atenuação de DEGs entre modelos |
| `auditoria_deconvolucao_qc.csv` | correlação células × marcadores |
| `auditoria_tss_estagio_child.csv` | TSS, estágio e Child-Pugh |
| `benchmarks_analises_avancadas.md` | tempos, memória e concordância das análises avançadas |
| `compliance_analises_avancadas.md` | checklist das análises avançadas |
| `manifesto_renomeacao.md` | regras de renomeação dos arquivos |

## Como reproduzir

```bash
Rscript scripts/etapa_10_auditoria_dados_volcano.R
Rscript scripts/etapa_11_benchmarks.R
Rscript scripts/etapa_12_compliance.R
```

Os relatórios são regerados a partir dos dados e figuras atuais.
