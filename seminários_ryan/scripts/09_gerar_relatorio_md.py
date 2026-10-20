#!/usr/bin/env python3
"""
09_gerar_relatorio_md.py
------------------------
Gera uma versão Markdown única de toda a documentação em docs/.

Uso (a partir de seminários_ryan/):
    python scripts/09_gerar_relatorio_md.py

Não requer dependências externas. O relatório final é somente Markdown;
não há geração de HTML.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
OUT = ROOT / "relatorio_completo.md"

ordem = [
    "01_analise_exploratoria.md",
    "02_ComBat.md",
    "03_logFC_limma.md",
    "04_volcano_plots.md",
    "05_relatorio_final.md",
    "06_dicionario_dados.md",
    "07_reproducibilidade.md",
    "08_fenotipos_avancados.md",
    "09_auditorias_benchmarks.md",
    "10_grupos_alternativos.md",
    "11_grupos_separados.md",
]

partes = []
for nome in ordem:
    p = DOCS / nome
    if not p.exists():
        continue
    texto = p.read_text(encoding="utf-8")
    # Na raiz, as imagens ficam em resultados/ (e não ../resultados/)
    texto = texto.replace("](../resultados/", "](resultados/")
    texto = texto.replace("](../auditorias/", "](auditorias/")
    # Links relativos entre documentos passam a apontar para docs/
    texto = re.sub(
        r"\]\((?!https?://|#|docs/)([^)]+\.md)\)",
        lambda m: f"](docs/{m.group(1)})",
        texto,
    )
    partes.append(texto.rstrip())

sumario = "\n".join(
    f"- [{nome.removesuffix('.md')}](docs/{nome})" for nome in ordem
    if (DOCS / nome).exists()
)

cabecalho = f"""# Relatório consolidado — Seminário Ryan (TCGA-LIHC × GTEx)

> Documento gerado a partir de `docs/*.md`.
> Versão **somente Markdown**; não há HTML consolidado.
> Para regenerar: `python scripts/09_gerar_relatorio_md.py`.

## Sumário

{sumario}
"""

OUT.write_text(cabecalho + "\n\n---\n\n" + "\n\n---\n\n".join(partes) + "\n",
               encoding="utf-8")
print(f"Gerado: {OUT}")
