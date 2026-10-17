#!/usr/bin/env python3
"""
09_gerar_relatorio_html.py
--------------------------
Gera uma versão HTML única e navegável de toda a documentação em docs/.

Uso (a partir de seminários_ryan/):
    python scripts/09_gerar_relatorio_html.py

Requisito: pip install markdown
"""
from pathlib import Path
import markdown

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
OUT = ROOT / "relatorio_completo.html"

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
]

partes = []
for nome in ordem:
    p = DOCS / nome
    if not p.exists():
        continue
    texto = p.read_text(encoding="utf-8")
    # No HTML na raiz, as imagens ficam em resultados/ (e não ../resultados/)
    texto = texto.replace("](../resultados/", "](resultados/")
    html = markdown.markdown(
        texto,
        extensions=["tables", "fenced_code", "toc", "sane_lists"],
    )
    partes.append(f'<section class="doc" id="{p.stem}">{html}</section>')

corpo = "\n<hr>\n".join(partes)

template = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Seminário Ryan — TCGA-LIHC × GTEx</title>
<style>
  :root { --azul:#1f4e79; --vermelho:#c0392b; --cinza:#555; }
  body { font-family: system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif;
         max-width: 1100px; margin: 2rem auto; padding: 0 1.2rem; line-height: 1.55;
         color: #222; background: #fff; }
  h1 { color: var(--azul); border-bottom: 2px solid #ddd; padding-bottom: .3rem; }
  h2 { color: var(--azul); margin-top: 2.2rem; border-bottom: 1px solid #eee; padding-bottom: .2rem; }
  h3 { color: #2c3e50; margin-top: 1.6rem; }
  a { color: var(--azul); }
  img { max-width: 100%; height: auto; display: block; margin: 1rem auto;
        border: 1px solid #eee; border-radius: 6px; }
  table { border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: .95rem; }
  th, td { border: 1px solid #ddd; padding: .45rem .6rem; vertical-align: top; }
  th { background: #f4f7fa; }
  tr:nth-child(even) { background: #fafafa; }
  code { background: #f4f4f4; padding: .1rem .3rem; border-radius: 3px; }
  pre { background: #f7f7f7; padding: .8rem; overflow-x: auto; border-radius: 6px; }
  blockquote { border-left: 4px solid #d0d7de; margin: 1rem 0; padding: .2rem 1rem;
               color: var(--cinza); background: #fafbfc; }
  hr { border: none; border-top: 2px solid #eee; margin: 2.5rem 0; }
  .doc { margin-bottom: 1rem; }
  .topo { background: #f4f7fa; border: 1px solid #dde5ee; border-radius: 8px;
          padding: 1rem 1.2rem; margin-bottom: 2rem; }
</style>
</head>
<body>
<div class="topo">
  <strong>Seminário Ryan — TCGA-LIHC × GTEx (fígado)</strong><br>
  Relatório consolidado gerado a partir de <code>docs/*.md</code>.
  Versão HTML para leitura; a fonte editável permanece em Markdown.
</div>
__CORPO__
</body>
</html>
"""

OUT.write_text(template.replace("__CORPO__", corpo), encoding="utf-8")
print(f"Gerado: {OUT}")
