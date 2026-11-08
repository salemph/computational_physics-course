from docx import Document
from docx.shared import Pt, RGBColor, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

d = Document("raw.docx")
BLUE = RGBColor(0x1C, 0x5C, 0xAB)

for sec in d.sections:
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.2)
    sec.top_margin = sec.bottom_margin = Cm(2.0)
    # page number footer
    p = sec.footer.paragraphs[0]
    p.alignment = 1
    run = p.add_run()
    for t, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        if t:
            e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), t); run._r.append(e)
        else:
            e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt; run._r.append(e)
    run.font.size = Pt(9)


def setfont(style, name, size=None, color=None, bold=None):
    f = style.font
    f.name = name
    rpr = style.element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts"); rpr.append(rf)
    for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rf.set(qn(a), name)
    for a in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
        if rf.get(qn(a)) is not None:
            del rf.attrib[qn(a)]
    if size: f.size = Pt(size)
    if color is not None: f.color.rgb = color
    if bold is not None: f.bold = bold


st0 = d.styles
class _S:
    def __getitem__(self, k):
        for x in st0:
            if x.name == k: return x
        raise KeyError(k)
    def __iter__(self): return iter(st0)
st = _S()
for nm in ("Normal", "Body Text", "First Paragraph", "Compact"):
    if nm in [s.name for s in st]:
        setfont(st[nm], "Calibri", 11)
        st[nm].paragraph_format.space_after = Pt(6)
        st[nm].paragraph_format.line_spacing = 1.12
setfont(st["Title"], "Calibri", 26, BLUE, True)
setfont(st["Subtitle"], "Calibri", 13, RGBColor(0x55, 0x55, 0x55), False)
for lvl, size in ((1, 18), (2, 14), (3, 12)):
    s = st[f"Heading {lvl}"]
    setfont(s, "Calibri", size, BLUE, True)
    s.paragraph_format.space_before = Pt(16 if lvl == 1 else 12)
    s.paragraph_format.space_after = Pt(6)
    s.paragraph_format.keep_with_next = True
    if lvl == 1:
        s.paragraph_format.page_break_before = True
for nm in ("Source Code",):
    if nm in [s.name for s in st]:
        s = st[nm]
        s.paragraph_format.space_after = Pt(0)
        s.paragraph_format.line_spacing = 1.0
        pPr = s.element.get_or_add_pPr()
        shd = OxmlElement("w:shd"); shd.set(qn("w:val"), "clear"); shd.set(qn("w:fill"), "F4F6F9")
        pPr.append(shd)
for s in st:
    if s.name in ("Verbatim Char",) or s.name.endswith("Tok"):
        setfont(s, "Consolas", 7.5)
for nm in ("Image Caption", "Table Caption", "Caption"):
    if nm in [s.name for s in st]:
        setfont(st[nm], "Calibri", 9.5, RGBColor(0x44, 0x44, 0x44))
        st[nm].font.italic = True
        st[nm].paragraph_format.space_after = Pt(10)

# code blocks: make sure every run in Source Code paragraphs is 7.5pt Consolas
for p in d.paragraphs:
    if p.style.name == "Source Code":
        for r in p.runs:
            r.font.size = Pt(7.5); r.font.name = "Consolas"
    if p.style.name in ("Captioned Figure", "Figure"):
        p.alignment = 1
    # centre images
    if p._p.xpath(".//pic:pic"):
        p.alignment = 1

# tables: borders, header shading, smaller font
for t in d.tables:
    tblPr = t._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        e = OxmlElement(f"w:{edge}"); e.set(qn("w:val"), "single"); e.set(qn("w:sz"), "4")
        e.set(qn("w:color"), "B7C4D6"); borders.append(e)
    tblPr.append(borders)
    tblW = tblPr.find(qn("w:tblW"))
    if tblW is None:
        tblW = OxmlElement("w:tblW"); tblPr.append(tblW)
    tblW.set(qn("w:w"), "5000"); tblW.set(qn("w:type"), "pct")
    for row_i, row in enumerate(t.rows):
        trPr = row._tr.get_or_add_trPr()
        cs = OxmlElement("w:cantSplit"); trPr.append(cs)
        if row_i == 0:
            th = OxmlElement("w:tblHeader"); trPr.append(th)
    for i, row in enumerate(t.rows):
        for cell in row.cells:
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(1)
                for r in p.runs:
                    r.font.size = Pt(8.5)
                    if i == 0: r.font.bold = True
            if i == 0:
                tcPr = cell._tc.get_or_add_tcPr()
                shd = OxmlElement("w:shd"); shd.set(qn("w:val"), "clear"); shd.set(qn("w:fill"), "E3ECF8")
                tcPr.append(shd)

settings = d.settings.element
uf = OxmlElement("w:updateFields"); uf.set(qn("w:val"), "true"); settings.append(uf)
d.save("Numerov_Schrodinger_Course.docx")
print("saved")
