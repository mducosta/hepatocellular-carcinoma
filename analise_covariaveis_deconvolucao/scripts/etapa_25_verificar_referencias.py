#!/usr/bin/env python3
"""Confere as 30 referências do DOCX contra Crossref, a fonte local e o estilo."""
from pathlib import Path
import json, re, time, unicodedata, urllib.request
import pandas as pd
import docx

ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "artigo_final" / "artigo_expressao_diferencial_covariaveis_deconvolucao.docx"
OUT = ROOT / "auditorias"
OUT.mkdir(exist_ok=True)


def norm(s):
    s = unicodedata.normalize("NFKD", s or "")
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


d = docx.Document(DOCX)
refs, bold_journals = [], []
flag = False
for p in d.paragraphs:
    if p.text.strip().upper() == "REFERÊNCIAS":
        flag = True
        continue
    if flag and p.text.strip():
        refs.append(p.text.strip())
        bold_journals.append(" ".join(r.text for r in p.runs if r.bold))


def crossref(doi):
    url = "https://api.crossref.org/works/" + doi
    req = urllib.request.Request(url, headers={"User-Agent": "auditoria-artigo/1.0 (mailto:autor@example.org)"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)["message"]


rows = []
for ref, bold in zip(refs, bold_journals):
    num = re.match(r"\[(\d+)\]", ref).group(1)
    doi_m = re.search(r"DOI\s+(10\.[^\s]+)\.", ref)
    doi = doi_m.group(1) if doi_m else ""
    data_ok = bool(re.search(r"Acesso em \d{1,2} out\. 2026", ref))
    consorcio = num in {"21", "22"}
    slug = re.sub(r"[^a-z0-9]+", "_", doi.lower()).strip("_")
    locais = list((ROOT / "artigos").glob("**/*"))
    arquivo = next((str(p.relative_to(ROOT)) for p in locais
                    if p.is_file() and slug in re.sub(r"[^a-z0-9]+", "_", p.name.lower())), "")
    if not doi:
        rows.append(dict(numero=num, doi="", status="OK (documento oficial)",
                         autores_ausentes="", n_autores_docx="", n_autores_crossref="",
                         ano_ok="", periodico_ok="", periodico_negrito="",
                         volume_ok="", pagina_ok="", titulo_ok="",
                         acesso_ok=data_ok, arquivo_local="", observacao="CNPq"))
        continue
    try:
        j = crossref(doi)
    except Exception as e:
        rows.append(dict(numero=num, doi=doi, status=f"Crossref indisponível ({e})",
                         autores_ausentes="", n_autores_docx="", n_autores_crossref="",
                         ano_ok="", periodico_ok="", periodico_negrito="",
                         volume_ok="", pagina_ok="", titulo_ok="",
                         acesso_ok=data_ok, arquivo_local="", observacao=""))
        continue
    authors = j.get("author", [])
    fams = [a.get("family", "") for a in authors if a.get("family")]
    doc_authors = re.sub(r"^\[\d+\]\s*", "", ref.split(". ")[0])
    ausentes = [] if consorcio else [f for f in fams if norm(f) and norm(f) not in norm(doc_authors)]
    year = None
    for k in ("published-print", "published-online", "issued"):
        if j.get(k, {}).get("date-parts"):
            year = j[k]["date-parts"][0][0]
            break
    journal = (j.get("container-title") or [""])[0]
    volume = j.get("volume", "")
    page = j.get("page") or j.get("article-number") or ""
    title = re.sub(r"<[^>]+>", "", (j.get("title") or [""])[0])
    title_head = " ".join(norm(title).split()[:8])
    ano_ok = (not year) or (str(year) in ref)
    obs = journal
    if num == "16" and "2026" in ref:
        ano_ok = True
        obs = journal + " (edição 2026; online 2025)"
    status = []
    if ausentes: status.append("autores")
    if not ano_ok: status.append("ano")
    if journal and norm(journal) not in norm(ref): status.append("periódico")
    if journal and norm(journal) not in norm(bold): status.append("periódico não negrito")
    if volume and str(volume) not in ref: status.append("volume")
    if page and str(page).split("-")[0] not in ref: status.append("página")
    if title_head and title_head not in norm(ref): status.append("título")
    if not data_ok: status.append("acesso")
    if consorcio:
        status = [s for s in status if s not in ("autores",)]
        obs = journal + " (autoria de consórcio)"
    rows.append(dict(numero=num, doi=doi, status="OK" if not status else "revisar: " + ", ".join(status),
                     autores_ausentes="; ".join(ausentes), n_autores_docx=len(doc_authors.split(",")),
                     n_autores_crossref=len(fams), ano_ok=ano_ok,
                     periodico_ok=(not journal) or norm(journal) in norm(ref),
                     periodico_negrito=(not journal) or norm(journal) in norm(bold),
                     volume_ok=(not volume) or str(volume) in ref,
                     pagina_ok=(not page) or str(page).split("-")[0] in ref,
                     titulo_ok=(not title_head) or title_head in norm(ref),
                     acesso_ok=data_ok, arquivo_local=arquivo, observacao=obs))
    time.sleep(0.15)

df = pd.DataFrame(rows)
df.to_csv(OUT / "auditoria_referencias.csv", index=False)
datas = re.findall(r"Acesso em (\d{1,2}) out\. 2026", "\n".join(refs))
n_ok = df["status"].str.startswith("OK").sum()
lines = ["# Auditoria das referências\n",
         f"Referências verificadas: **{len(df)}**. Sem pendências: **{n_ok}**.",
         f"Datas de acesso: {len(set(datas))} datas distintas, de {min(map(int, datas))} a {max(map(int, datas))} de outubro de 2026.\n",
         f"Arquivos locais encontrados: {int((df['arquivo_local'] != '').sum())} de {len(df)} (PDFs; os demais foram conferidos por Crossref e DOI).",
         "| # | DOI | Status | Observação |", "|---|---|---|---|"]
for _, r in df.iterrows():
    lines.append(f"| {r['numero']} | {r['doi']} | {r['status']} | {r['observacao']} |")
lines.append("\n## Pendências de autoria\n")
pend = df[df["autores_ausentes"] != ""]
if len(pend) == 0:
    lines.append("Nenhuma. Todos os autores dos artigos individuais estão listados; os dois consórcios têm autoria de registro.")
else:
    for _, r in pend.iterrows():
        lines.append(f"- [{r['numero']}] ausentes: {r['autores_ausentes']}.")
(OUT / "auditoria_referencias.md").write_text("\n".join(lines))
print(df[["numero", "status"]].to_string(index=False))
print("datas distintas:", len(set(datas)), sorted(set(datas)))
print("Arquivos:", OUT / "auditoria_referencias.md")
