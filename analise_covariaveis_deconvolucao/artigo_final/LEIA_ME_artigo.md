# Artigo final

**Título:** Expressão diferencial no hepatocarcinoma e o que permanece após
ajustar por clínica, lote e composição celular

**Autor:** Ryan Santos

**Espécie:** *Homo sapiens* (tecido hepático).

**Dados:** TCGA-LIHC e GTEx via UCSC Xena; expressão Toil v2016-09-03,
fenótipos GTEx v8; acesso em 17/10/2026.

**Versões:** destacadas na seção inicial do DOCX (bases de dados, R 4.6.1,
Bioconductor 3.23, pacotes e Python 3.14.5).

**Artigo:** 26/10/2026.

## Arquivos

- `artigo_expressao_diferencial_covariaveis_deconvolucao.docx` — versão editável.
- `figuras/` — figuras em PNG, 300 dpi, usadas no artigo.

## Formatação

A4; Times New Roman 12; texto justificado; títulos em negrito; tabelas com
grade; legendas de figuras e tabelas com fonte; referências numeradas com
todos os autores, periódico em negrito e datas de acesso distintas. Sem resumo,
introdução, cabeçalho, rodapé ou paginação.

As versões das bases de dados, do ambiente computacional e dos pacotes estão
destacadas na seção inicial do artigo e detalhadas em `../docs/`.

## Observações

- O artigo usa 35 referências; todas as citações do corpo têm referência
  correspondente. As datas de acesso vão de 15 a 28 de outubro de 2026.
- O texto não usa marcadores de lista nem dois pontos em títulos, subtítulos ou
  frases do corpo. Termos estrangeiros estão em itálico.
- A discussão inclui a subseção "O que este estudo acrescenta à literatura".
- A Figura 5 foi refeita. O painel F, que mostra o modelo de estágio restrito
  aos mesmos casos completos, agora exibe o limiar de FDR, a contagem zero de
  DEGs e a anotação do FDR mínimo (0,64), evitando a impressão de painel vazio.
- Todas as figuras passaram por verificação automática de sobreposição de
  rótulos. A Figura 5 usa posicionamento com linhas-guia, e as Figuras 8 e 9
  tiveram rótulos encurtados, rotação e mais espaço entre painéis. As tabelas 4
  e 5 foram reajustadas para caber nos 17 cm úteis da página.
- O artigo de Kočar et al. (2026), que justifica o dimorfismo sexual em
  hepatócitos, está entre as referências citadas.
- Os PDFs/XMLs/HTMLs das referências são mantidos apenas localmente; a pasta
  `../artigos/` não é versionada.

## Análises de sensibilidade e validação

O artigo inclui oito tabelas suplementares. A Tabela S1 traz a sensibilidade do
ajuste por composição, a S2 a ponderação de probabilidade inversa, a S3 a
validação em quatro coortes independentes (GSE14520, GSE36376, GSE76427 e
GSE62232), a S4 a validação interna do grau, a S5 os modelos mistos com TSS, a
S6 a comparação com contagens brutas (DESeq2 e edgeR), a S7 os 48 modelos e a S8
os 27 genes de grau preservados pelo escore curado.

As auditorias finais estão em `../auditorias/auditoria_numerica.md`,
`../auditorias/auditoria_referencias.md`,
`../auditorias/auditoria_child_pugh_gdc.md` e
`../auditorias/auditoria_gdc_outros.md`. A auditoria numérica cobre 137
verificações, e as 35 referências foram conferidas contra Crossref e os PDFs
locais.
