# AUDITORIA DO RESUMO EXPANDIDO — CRUZAMENTO COM SCRIPT E DADOS

**Arquivo:** resumo_expandido_hepatocarcinoma_maria_costa.docx  
**Data da auditoria:** 2026-06-26  
**Método:** Comparação de cada afirmação com (a) código R corrigido, (b) arquivo de dados disponível, (c) genes da via KEGG hsa05417, (d) literatura citada

---

## CLASSIFICAÇÃO DAS AFIRMAÇÕES

Legenda:
- ✅ CORRETA
- ⚠️ PARCIALMENTE CORRETA
- ❌ NÃO SUPORTADA (dados indisponíveis)
- 🔴 FALSA / ERRO CRÍTICO
- 🟡 EXAGERADA
- 🔵 ESPECULATIVA

---

## SEÇÃO: TÍTULO

> "Reprogramação Transcriptômica da Via de Lipídios e Aterosclerose no Hepatocarcinoma"

| Classificação | Justificativa |
|---------------|---------------|
| ⚠️ PARCIALMENTE CORRETA | O título sugere que o estudo DEMONSTRA reprogramação. Seria mais preciso: "Análise Exploratória da Expressão de Genes da Via..." ou "Perfil Transcriptômico..." |

---

## SEÇÃO: RESUMO

### Afirmação 1
> "Este estudo transcriptômico exploratório utilizou dados públicos do TCGA e GTEx para identificar genes diferencialmente expressos (DEGs) na via de Lipídios e Aterosclerose em amostras de CHC."

| Classificação | Justificativa |
|---------------|---------------|
| ⚠️ PARCIALMENTE CORRETA | O desenho do estudo está correto, mas os dados de expressão gênica NÃO estão disponíveis no arquivo atual. A identificação dos DEGs não pôde ser executada. |

### Afirmação 2
> "As análises foram realizadas no ambiente R com o pacote limma, adotando FDR < 0,05 e logFC > 1 como critérios de significância."

| Classificação | Justificativa |
|---------------|---------------|
| ✅ CORRETA | O script usa limma com FDR < 0.05 e |logFC| > 1. Metodologia estatística adequada. |

### Afirmação 3
> "Entre os DEGs identificados, destacam-se CALM2, associada à proliferação tumoral e pior prognóstico; CAMK4, regulador da proliferação de células hepáticas parenquimatosas; MMP9, promotora da invasão tumoral e da resistência a terapias sistêmicas; e BAX, cujo comprometimento favorece a evasão apoptótica e a progressão neoplásica."

| Gene | Classificação | Justificativa |
|------|---------------|---------------|
| CALM2 | ❌ NÃO SUPORTADA | Gene está na via KEGG hsa05417, mas expressão diferencial NÃO pode ser confirmada (sem dados de expressão). Afirmação de "pior prognóstico" é EXAGERADA sem validação clínica. |
| CAMK4 | 🔴 ERRO CRÍTICO | **CAMK4 NÃO pertence à via KEGG hsa05417 (Lipid and Atherosclerosis)!** A via contém CAMK2A, CAMK2B, CAMK2D, CAMK2G, mas NÃO CAMK4. O resumo discute um gene que NÃO está na via analisada. |
| MMP9 | ❌ NÃO SUPORTADA | Expressão diferencial não pode ser confirmada. Afirmações sobre "resistência a terapias sistêmicas" são ESPECULATIVAS. |
| BAX | ❌ NÃO SUPORTADA | Expressão diferencial não pode ser confirmada. Afirmações sobre "potencial alvo terapêutico" são ESPECULATIVAS. |

### Afirmação 4
> "Clinicamente, o CHC exige vigilância periódica, diagnóstico precoce e tratamento individualizado, com destaque para imunoterápicos e inibidores de tirosina quinase."

| Classificação | Justificativa |
|---------------|---------------|
| ⚠️ PARCIALMENTE CORRETA | Informação clínica genérica correta, mas desconectada dos resultados do estudo (que são transcriptômicos, não clínicos). Menciona tratamentos sem relação com os genes analisados. |

### Afirmação 5
> "A enfermagem oncológica desempenha papel essencial na coordenação do cuidado e na detecção precoce da doença."

| Classificação | Justificativa |
|---------------|---------------|
| ⚠️ PARCIALMENTE CORRETA | Papel da enfermagem é relevante, mas a afirmação é genérica e não conecta com os achados transcriptômicos do estudo. |

---

## SEÇÃO: 1. INTRODUÇÃO

### Afirmação 6
> "O carcinoma hepatocelular (CHC) constitui uma das neoplasias de maior relevância clínica e epidemiológica global, sendo a terceira principal causa de mortalidade associada ao câncer no mundo[1], com prognóstico desfavorável especialmente nos estágios avançados[2]."

| Classificação | Justificativa |
|---------------|---------------|
| ✅ CORRETA | Informação epidemiológica precisa. Referências [1] e [2] são compatíveis. |

### Afirmação 7
> "Quinbin et al., por meio de abordagens transcriptômicas, metabolômicas e lipidômicas, identificaram alterações em vias metabólicas [...] evidenciando a reprogramação metabólica celular como característica fundamental da oncogênese hepática[3]."

| Classificação | Justificativa |
|---------------|---------------|
| ⚠️ PARCIALMENTE CORRETA | O nome correto do primeiro autor é **Liu Q** (não "Quinbin"). A referência [3] é "LIU, Q.; ZHANG, X.; QI, J. et al. Comprehensive profiling of lipid metabolic reprogramming..." — o resumo cita "Quinbin et al." quando deveria ser "Liu et al." |

### Afirmação 8
> "a enfermagem desempenha papel essencial na tradução desse conhecimento à prática clínica, por meio da vigilância contínua, promoção da adesão ao rastreamento e identificação precoce de alterações sugestivas de progressão da doença em pacientes de risco, como portadores de cirrose hepática[5]."

| Classificação | Justificativa |
|---------------|---------------|
| 🟡 EXAGERADA | A referência [5] (Oikonomou et al., 2025) trata de intervenção educacional em cirrose, não especificamente de enfermagem em genética/genômica ou interpretação de dados transcriptômicos. A conexão com o estudo é tênue. |

---

## SEÇÃO: 2. MATERIAIS E MÉTODOS

### Afirmação 9
> "Foram utilizados dados públicos de expressão gênica e informações clínicas de amostras de CHC obtidos do The Cancer Genome Atlas (TCGA) e amostras de tecido hepático normal provenientes do Genotype-Tissue Expression Project (GTEx), ambos acessados por meio da plataforma UCSC Xena."

| Classificação | Justificativa |
|---------------|---------------|
| ⚠️ PARCIALMENTE CORRETA | A descrição da fonte de dados está correta, mas o arquivo atual contém APENAS metadados, sem expressão gênica. |

### Afirmação 10
> "Os valores de expressão gênica previamente normalizados foram analisados em escala log2."

| Classificação | Justificativa |
|---------------|---------------|
| ❌ NÃO SUPORTADA | Sem dados de expressão gênica, não é possível confirmar a escala. O script detecta automaticamente se usa log2 (limma direto) ou contagem (voom). |

### Afirmação 11
> "A identificação dos genes diferencialmente expressos (DEGs) foi realizada utilizando o pacote limma por meio de modelagem linear e correção para múltiplos testes pelo método de Benjamini-Hochberg (FDR/adj.P.Val). Foram considerados diferencialmente expressos os genes que apresentaram log Fold Change > 1 (logFC) e taxa de falsos descobrimentos (FDR) <0,05."

| Classificação | Justificativa |
|---------------|---------------|
| ✅ CORRETA | Metodologia descrita corresponde exatamente ao script corrigido. |

### Afirmação 12
> "os genes significativos foram submetidos à construção de rede de interação proteína-proteína (PPI) na plataforma STRING, considerando escore de confiança ≥ 700, bem como da construção de um volcano plot com vistas a demonstrar os genes super/subexpressos, conforme seu logFC e significância."

| Classificação | Justificativa |
|---------------|---------------|
| ✅ CORRETA | O script novo usa STRINGdb v12.0 com score ≥ 700 e network_type = "physical". O volcano plot foca nos genes da via LA. |

---

## SEÇÃO: 3. RESULTADOS E DISCUSSÃO

### Afirmação 13 (CRÍTICA)
> "Os resultados da análise transcriptômica, evidenciados no volcano plot, identificaram genes relacionados à via de lipídios e aterosclerose, com destaque para as calmodulinas (CALML5, CALM2), CAMK4, BAX e MMP9."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | **Resultados NÃO podem ser apresentados como fato**, pois a análise NÃO foi executada (sem dados de expressão gênica). Além disso, **CAMK4 NÃO está na via KEGG hsa05417**. |

### Afirmação 14
> "Entre os DEGs, CALML5 apresentou maior relevância estatística; contudo, a literatura não estabelece associação direta entre sua expressão e a patogênese do hepatocarcinoma."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | **Dados inventados.** A análise não foi executada. Não é possível afirmar que CALML5 "apresentou maior relevância estatística". |

### Afirmação 15
> "Diferentemente, CALM2, membro da família das calmodulinas e principal reguladora da sinalização dependente de Ca²⁺, modula a proliferação e migração de células tumorais, associando-se a pior prognóstico clínico[6]. Seu silenciamento por siRNA em linhagens de carcinoma hepatocelular inibiu a proliferação e induziu apoptose[7]."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | **Afirmações de prognóstico sem validação.** O estudo é exploratório com dados secundários. Não pode estabelecer associação prognóstica. As referências [6] e [7] são artigos legítimos, mas a conexão "CALM2 → pior prognóstico" extrapola os resultados possíveis de uma análise transcriptômica exploratória. |

### Afirmação 16
> "No que se refere ao CAMK4, gene-alvo a jusante de CAMKK2, este integra uma via de sinalização responsável pela regulação da proliferação de células parenquimatosas hepáticas[8,9]."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | **CAMK4 NÃO pertence à via KEGG hsa05417.** Discutir CAMK4 como resultado da análise da via LA é cientificamente incorreto. Se o gene não está na via, não deveria ser mencionado como resultado do estudo. |

### Afirmação 17
> "No que concerne à MMP9, sua superexpressão promove a degradação da matriz extracelular, incluindo a membrana basal, favorecendo a invasão tumoral e a remodelação do microambiente imunossupressor, fatores associados à resistência a terapias sistêmicas[10,11]."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | **Afirmações de mecanismo molecular confirmado e resistência terapêutica** sem suporte de análise funcional. O estudo é transcriptômico exploratório — não demonstra mecanismo, não demonstra resistência terapêutica. |

### Afirmação 18
> "Por fim, BAX, membro pró-apoptótico da família BCL-2, desempenha papel fundamental na permeabilização da membrana mitocondrial externa (MOMP); sua desregulação compromete a ativação da cascata apoptótica, conferindo vantagem proliferativa às células neoplásicas e consolidando-o como potencial alvo terapêutico[12]."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | **"Potencial alvo terapêutico"** é afirmação que requer validação funcional e clínica. O estudo exploratório não pode fazer este tipo de afirmação. |

### Afirmação 19
> "Nesse sentido, torna-se necessário compreender o CHC para além de sua biologia molecular, abrangendo a perspectiva clínica. Por ser assintomático nos estágios iniciais, programas de vigilância são essenciais [...]"

| Classificação | Justificativa |
|---------------|---------------|
| ⚠️ PARCIALMENTE CORRETA | A transição para conteúdo clínico é abrupta e desconectada dos resultados transcriptômicos. O parágrafo inteiro sobre vigilância, APASL, diagnóstico por imagem, etc. é clinicamente correto mas NÃO é resultado do estudo. |

### Afirmação 20
> "Estudo conduzido por Nazareth et al. avaliou a eficácia de clínica de vigilância do CHC gerenciada por enfermeiros, demonstrando melhora na adesão, detecção precoce e otimização da carga assistencial[13]."

| Classificação | Justificativa |
|---------------|---------------|
| ✅ CORRETA | Referência compatível: NAZARETH, S. et al. Nurse-led hepatocellular carcinoma surveillance clinic. Int J Nurs Pract, 2016. |

### Afirmação 21
> "Tal abordagem articula-se com o Modelo de Adaptação de Callista Roy, que preconiza o acompanhamento contínuo do paciente crônico a partir dos modos adaptativos fisiológico, autoconceito, função de papel e interdependência[16]."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | A referência [16] é "CRAIG, Amanda J. et al. Tumour evolution in hepatocellular carcinoma. Nature Reviews Gastroenterology & Hepatology, v. 17, p. 139–152, 2020." Este artigo é sobre evolução tumoral do CHC, **NÃO sobre o Modelo de Adaptação de Callista Roy**. A citação é INCOMPATÍVEL com o conteúdo. |

---

## SEÇÃO: 4. CONCLUSÕES

### Afirmação 22
> "Em conjunto, os resultados encontrados evidenciam que a compreensão da biologia do hepatocarcinoma deve ser integrada ao acompanhamento clínico contínuo, incluindo a relevância da expressão gênica na compreensão dos mecanismos de progressão tumoral e na tomada de decisão clínica, no qual a enfermagem desempenha papel central na coordenação do cuidado e na identificação precoce de alterações clínicas associadas à evolução da doença, contribuindo para a qualificação da assistência e para melhores desfechos terapêuticos."

| Classificação | Justificativa |
|---------------|---------------|
| 🔴 ERRO CRÍTICO | **"Tomada de decisão clínica"** baseada em análise transcriptômica exploratória é afirmação não suportada. **"Melhores desfechos terapêuticos"** é alegação de efetividade sem validação. A conclusão extrapola massivamente os limites de um estudo exploratório com dados secundários. |

---

## RESUMO DOS ERROS CRÍTICOS

| # | Erro | Severidade |
|---|------|------------|
| 1 | **CAMK4 não pertence à via KEGG hsa05417** | CRÍTICO |
| 2 | Resultados apresentados como fato sem análise executada | CRÍTICO |
| 3 | "CALML5 apresentou maior relevância estatística" — dado inventado | CRÍTICO |
| 4 | Alegações de prognóstico (CALM2 → pior prognóstico) | CRÍTICO |
| 5 | Alegações de mecanismo molecular confirmado (MMP9) | CRÍTICO |
| 6 | Alegações de alvo terapêutico (BAX) | CRÍTICO |
| 7 | Alegações de resistência terapêutica (MMP9) | CRÍTICO |
| 8 | Referência [16] incompatível (Craig ≠ Callista Roy) | CRÍTICO |
| 9 | Citação "Quinbin et al." — nome incorreto (deveria ser Liu et al.) | MODERADO |
| 10 | Conclusão extrapola para "decisão clínica" e "desfechos terapêuticos" | CRÍTICO |

---

## VEREDITO DO RESUMO EXPANDIDO

| Critério | Nota |
|----------|------|
| Coerência científica | 2/10 🔴 |
| Precisão dos resultados | 1/10 🔴 |
| Referências | 4/10 🟡 |
| Conexão com enfermagem | 5/10 ⚠️ |
| Prudência científica | 1/10 🔴 |
| Adequação ao template | 10/10 ✅ |
| **Status geral** | **⛔ NÃO APROVADO — REQUER CORREÇÃO SUBSTANCIAL** |
