#!/usr/bin/env python3
# =============================================================================
# M36_build_word_ajhg_20260816.py — Assemble the AJHG-format submission Word
# document from the verified manuscript (manuscript.md), references, and the
# English-label figures.
#   Output: docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx
# Formatting per AJHG author guidelines (2026-08-16 revision):
#   - No separate title page (author instruction): title + author + affiliation + running title
#   - Times New Roman, 12 pt, double-spaced main text
#   - Numeric superscript citations ([@key] -> superscript number)
#   - Abstract <=200 words
#   - IMRaD sections, then Data and Code Availability, Web Resources,
#     Declaration of Interests, Funding
#   - References in NLM style (numbered by first citation order)
#   - Embedded figures (9 + S1), each with its full legend beneath (AJHG style)
#   - One three-line table (Table 1) with caption and note
# =============================================================================
import json, os, re

from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
DOCS = f"{BASE}/docs/manuscript"
FIGD = f"{BASE}/results/figures"
OUTP = f"{DOCS}/AJHG_submission_Qiushuo_Geng_20260816.docx"

REFS = json.load(open(f"{DOCS}/refs.json"))

FIG_FILES = {
    1: f"{FIGD}/20260816_Fig1_design_genome.png",
    2: f"{FIGD}/20260816_Fig2_yield_funnel.png",
    3: f"{FIGD}/20260816_Fig3_calibration.png",
    4: f"{FIGD}/20260816_Fig4_regional.png",
    5: f"{FIGD}/20260816_Fig5_candidates.png",
}
S1_FILE = f"{FIGD}/20260816_FigS1_susie.png"

# combined inline tokenizer: **bold**, *italic*, [@key]
INLINE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|\[@[a-z0-9]+\])")

def set_style_base(doc):
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(12)
    style.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    rpr = style.element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts"); rpr.append(rfonts)
    rfonts.set(qn("w:ascii"), "Times New Roman")
    rfonts.set(qn("w:hAnsi"), "Times New Roman")
    rfonts.set(qn("w:eastAsia"), "Times New Roman")

def add_run(p, text, bold=False, italic=False, superscript=False, font_size=12):
    r = p.add_run(text)
    r.font.name = "Times New Roman"
    r.font.size = Pt(font_size)
    r.font.bold = bold
    r.font.italic = italic
    if superscript:
        r.font.superscript = True
    return r

def set_double(p, after=0):
    pf = p.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    pf.space_after = Pt(after)

def add_inline(p, text, font_size=12, base_bold=False):
    """Add text to an existing paragraph, handling **bold**, *italic*, [@key]."""
    pos = 0
    for m in INLINE.finditer(text):
        if m.start() > pos:
            add_run(p, text[pos:m.start()], bold=base_bold, font_size=font_size)
        tok = m.group(1)
        if tok.startswith("[@"):
            key = tok[2:-1]
            num = REFS.get(key, "?")
            if num == "?": print(f"  !! missing citation key: {key}")
            add_run(p, str(num), bold=base_bold, superscript=True, font_size=font_size)
        elif tok.startswith("**"):
            add_run(p, tok[2:-2], bold=True, font_size=font_size)
        elif tok.startswith("*"):
            add_run(p, tok[1:-1], italic=True, font_size=font_size)
        pos = m.end()
    if pos < len(text):
        add_run(p, text[pos:], bold=base_bold, font_size=font_size)

def add_body_para(doc, text, first_indent=None):
    p = doc.add_paragraph(); set_double(p)
    if first_indent: p.paragraph_format.first_line_indent = Inches(first_indent)
    add_inline(p, text)

def add_heading(doc, text, level=1):
    p = doc.add_paragraph(); set_double(p)
    if level == 1:
        add_run(p, text, bold=True, font_size=13)
    else:
        add_run(p, text, bold=True, font_size=12)
    p.paragraph_format.keep_with_next = True
    return p

def set_cell_borders(cell, top=None, bottom=None):
    tcPr = cell._tc.get_or_add_tcPr()
    borders = tcPr.find(qn("w:tcBorders"))
    if borders is None:
        borders = OxmlElement("w:tcBorders"); tcPr.append(borders)
    def set_edge(name, sz, val="single"):
        el = borders.find(qn(f"w:{name}"))
        if el is None:
            el = OxmlElement(f"w:{name}"); borders.append(el)
        el.set(qn("w:val"), val)
        el.set(qn("w:sz"), str(sz))
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), "000000")
    if top is not None: set_edge("top", top)
    if bottom is not None: set_edge("bottom", bottom)

def add_three_line_table(doc, header, rows, font_size=10):
    t = doc.add_table(rows=1 + len(rows), cols=len(header))
    t.alignment = WD_TABLE_ALIGNMENT.LEFT
    t.autofit = False
    for j, h in enumerate(header):
        cell = t.cell(0, j)
        cell.text = ""
        p = cell.paragraphs[0]
        p.paragraph_format.space_after = Pt(2)
        add_inline(p, h, font_size=font_size, base_bold=True)
        set_cell_borders(cell, top=12, bottom=6)
    for i, row in enumerate(rows, start=1):
        last = (i == len(rows))
        for j, v in enumerate(row):
            cell = t.cell(i, j)
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(2)
            add_inline(p, v, font_size=font_size)
            set_cell_borders(cell, bottom=(12 if last else 0))
    return t

def add_caption(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(14)
    p.paragraph_format.space_after = Pt(4)
    set_double(p)
    add_inline(p, text, font_size=11, base_bold=True)

def add_note(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(10)
    set_double(p)
    add_inline(p, text, font_size=10)

def add_figure_page(doc, fig_num, img_path, legend_text):
    doc.add_page_break()
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(); run.add_picture(img_path, width=Inches(6.3))
    cap = doc.add_paragraph(); set_double(cap)
    add_inline(cap, f"**Figure {fig_num}.** " + legend_text, font_size=11)

def parse_refs():
    refs = {}
    for line in open(f"{DOCS}/references.md", encoding="utf-8"):
        m = re.match(r"^(\d+)\.\s+(.*)", line.strip())
        if m: refs[int(m.group(1))] = m.group(2).strip()
    return refs

def split_sections(md_text):
    """Split manuscript.md into named sections; return list of (section_name, lines).
    Section boundaries are lines starting with '## '."""
    lines = md_text.split("\n")
    sections = []
    cur_name = None; cur = []
    for line in lines:
        m = re.match(r"^##\s+(.*)", line)
        if m:
            if cur_name is not None:
                sections.append((cur_name, cur))
            cur_name = m.group(1).strip(); cur = []
        else:
            if cur_name is not None:
                cur.append(line)
    if cur_name is not None:
        sections.append((cur_name, cur))
    return sections

def parse_paragraph_blocks(lines):
    """Yield ('h3'|'table'|'p'|'list'|'blank', payload) from section body lines."""
    blocks = []
    in_table = False; table_rows = []
    for line in lines:
        raw = line.rstrip()
        s = raw.strip()
        if not s:
            if in_table:
                blocks.append(("table", table_rows)); table_rows = []; in_table = False
            continue
        if s in ("---", "***", "___"):
            continue
        if s.startswith("|"):
            if not in_table: in_table = True; table_rows = []
            if not re.match(r"^\|[\s\-:|]+\|$", s):
                table_rows.append([c.strip() for c in s.strip("|").split("|")])
            continue
        if in_table:
            blocks.append(("table", table_rows)); table_rows = []; in_table = False
        if s.startswith("###"):
            blocks.append(("h3", s[3:].strip()))
            continue
        if s.startswith("- ") or s.startswith("* "):
            blocks.append(("list", s[2:].strip()))
            continue
        blocks.append(("p", s))
    if in_table:
        blocks.append(("table", table_rows))
    return blocks

def strip_figure_legend(text):
    """From a '**Figure N. <Title>.** <rest>' paragraph, return (num, '<Title>. <rest>').
    The 'Figure N.' prefix is added by add_figure_page, so it is NOT included here."""
    m = re.match(r"^\*\*Figure (\d+)\.\s*(.*?)\*\*(.*)$", text)
    if not m:
        return None, text
    num = int(m.group(1))
    title = m.group(2).strip().rstrip(".")
    rest = m.group(3).strip()
    body = f"{title}."
    if rest: body += " " + rest
    return num, body

def main():
    doc = Document()
    set_style_base(doc)
    for section in doc.sections:
        section.top_margin = Inches(1); section.bottom_margin = Inches(1)
        section.left_margin = Inches(1); section.right_margin = Inches(1)

    md_text = open(f"{DOCS}/manuscript.md", encoding="utf-8").read()
    sections = split_sections(md_text)
    by_name = {name: lines for name, lines in sections}

    # ---------- Title (no separate title page per author instruction) ----------
    title = "A transcriptome-wide cis-MR and colocalization atlas for type 2 diabetes, coronary artery disease, and fasting glucose: operating characteristics and candidate effector genes"
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; set_double(p)
    add_run(p, title, bold=True, font_size=14)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; set_double(p)
    add_run(p, "Qiushuo Geng", bold=True)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; set_double(p)
    add_run(p, "School of Medical Devices, Shenyang Pharmaceutical University, Benxi 117004, China")
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; set_double(p)
    add_run(p, "Running title: Transcriptome-wide cis-MR × coloc atlas for metabolic traits")
    doc.add_paragraph()

    # ---------- Main body (Abstract, numbered sections, data availability, web/declaration/funding) ----------
    BODY_ORDER = ["Abstract", "1. Introduction", "2. Methods", "3. Results",
                  "4. Discussion", "Data and Code Availability", "Web Resources",
                  "Declaration of Interests", "Funding"]
    for name in BODY_ORDER:
        if name not in by_name:
            print(f"  !! section not found: {name}")
            continue
        if name == "Abstract":
            add_heading(doc, "Abstract", 1)
        elif name.startswith(("1.", "2.", "3.", "4.")):
            add_heading(doc, name, 1)
        else:
            add_heading(doc, name, 1)
        for kind, payload in parse_paragraph_blocks(by_name[name]):
            if kind == "h3":
                add_heading(doc, payload, 2)
            elif kind == "p":
                add_body_para(doc, payload)
            elif kind == "list":
                p = doc.add_paragraph(); set_double(p)
                p.paragraph_format.left_indent = Inches(0.35)
                add_inline(p, payload)
            elif kind == "table":
                print("  !! table inside body section (unexpected):", payload[0][:4])

    # ---------- References ----------
    doc.add_page_break()
    add_heading(doc, "References", 1)
    refs_list = parse_refs()
    for num in sorted(refs_list):
        p = doc.add_paragraph(); set_double(p)
        p.paragraph_format.left_indent = Inches(0.4)
        p.paragraph_format.first_line_indent = Inches(-0.4)
        add_run(p, f"{num}. "); add_run(p, refs_list[num])

    # ---------- Tables (from 'Tables' section) ----------
    if "Tables" in by_name:
        doc.add_page_break()
        add_heading(doc, "Tables", 1)
        for kind, payload in parse_paragraph_blocks(by_name["Tables"]):
            if kind == "p":
                add_caption(doc, payload)
            elif kind == "table":
                add_three_line_table(doc, payload[0], payload[1:])
            elif kind == "list":
                add_note(doc, payload)

    # ---------- Figure Legends + embedded figures ----------
    figure_legends = {}
    if "Figure Legends" in by_name:
        for kind, payload in parse_paragraph_blocks(by_name["Figure Legends"]):
            if kind == "p":
                num, body = strip_figure_legend(payload)
                if num is not None:
                    figure_legends[num] = body
    for num in range(1, 6):
        legend = figure_legends.get(num, "")
        if not legend:
            print(f"  !! missing legend for Figure {num}")
        img = FIG_FILES[num]
        if os.path.exists(img):
            add_figure_page(doc, num, img, legend)
        else:
            print(f"  !! missing figure image: {img}")

    # Supplemental Fig S1 (resource snapshot)
    doc.add_page_break()
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(); run.add_picture(S1_FILE, width=Inches(6.3))
    cap = doc.add_paragraph(); set_double(cap)
    add_inline(cap, "**Supplemental Figure S1.** coloc.susie sensitivity of strong-colocalization calls. "
                    "Paired coloc.abf vs coloc.susie PP.H4 for the six adjudicated loci (RBM6×T2D, CNNM2×CAD, "
                    "PLAUR×CAD, CD101×T2D, RIC8A×CAD, LAMC1×CAD). SuSiE credible-set counts per side and "
                    "non-convergence markers (✗, all runs) are annotated. LAMC1×CAD (abf PP.H4 = 0.9139 → "
                    "SuSiE 0.0000) is highlighted; it is excluded from the candidate set on independent FDR and "
                    "multi-signal grounds. Because coloc.susie did not converge under external LD for any tested "
                    "locus, SuSiE posteriors are reported as exploratory only and are not used for primary inference.",
               font_size=11)

    # ---------- Supplemental Items ----------
    doc.add_page_break()
    add_heading(doc, "Supplemental Items", 1)
    if "Supplemental Items" in by_name:
        for kind, payload in parse_paragraph_blocks(by_name["Supplemental Items"]):
            if kind in ("p", "list"):
                p = doc.add_paragraph(); set_double(p)
                p.paragraph_format.left_indent = Inches(0.35)
                add_inline(p, payload)

    doc.save(OUTP)
    print("Saved:", OUTP)

if __name__ == "__main__":
    main()
