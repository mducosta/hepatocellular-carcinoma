#!/usr/bin/env python3
"""
etapa_09_gerar_relatorio_consolidado.py
------------------------
Gera uma versão Markdown única de toda a documentação em docs/.

Uso (a partir de seminários_ryan/):
    python scripts/etapa_09_gerar_relatorio_consolidado.py

Não requer dependências externas. O relatório final é somente Markdown;
não há geração de HTML.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
OUT = ROOT / "relatorio_consolidado_seminario.md"

ordem = [
    "analise_exploratoria_fenotipos.md",
    "correcao_lote_ComBat.md",
    "calculo_logFC_limma.md",
    "volcano_plots.md",
    "relatorio_final.md",
    "dicionario_dados.md",
    "reprodutibilidade.md",
    "fenotipos_clinicos_avancados.md",
    "auditorias_benchmarks.md",
    "grupos_alternativos.md",
    "analises_grupos_separados.md",
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
        r"\]\((?!https?://|#|docs/|auditorias/|resultados/)([^)]+\.md)\)",
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
> Para regenerar: `python scripts/etapa_09_gerar_relatorio_consolidado.py`.

## Sumário

{sumario}
"""

OUT.write_text(cabecalho + "\n\n---\n\n" + "\n\n---\n\n".join(partes) + "\n",
               encoding="utf-8")
print(f"Gerado: {OUT}")
