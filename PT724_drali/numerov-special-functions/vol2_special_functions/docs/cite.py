"""Numbered-citation processor for the course markdown.

Syntax in the text:  [@key]   [@key, locator]   [@k1; @k2, loc]   [@k1; prefix text @k2]
* numbers are assigned in order of first appearance in the whole document;
* every top-level section ("# ...") that cites something gets a
  "References for this chapter" list appended at its end;
* the placeholder {{BIBLIOGRAPHY}} is replaced by the complete numbered list,
  followed by "Further reading" (entries of the bibliography never cited).
"""
import re

CITE = re.compile(r"\[(@[^\[\]]*?|[^\[\]@]*?@[^\[\]]*?)\]")
PART = re.compile(r"^(.*?)@([A-Za-z0-9_]+)(?:,\s*(.*))?$")


def _parts(body):
    out = []
    for p in body.split(";"):
        m = PART.match(p.strip())
        if not m:
            raise ValueError(f"bad citation part: {p!r}")
        out.append((m.group(1).strip(), m.group(2), (m.group(3) or "").strip()))
    return out


def process(md, refs, chapter_heading="References for this chapter"):
    order = {}
    for m in CITE.finditer(md):
        for _, key, _ in _parts(m.group(1)):
            if key not in refs:
                raise KeyError(f"unknown reference key: {key}")
            order.setdefault(key, len(order) + 1)

    def render(m):
        items = []
        for pre, key, loc in _parts(m.group(1)):
            s = (pre + " " if pre else "") + str(order[key]) + (", " + loc if loc else "")
            items.append(s)
        return "[" + "; ".join(items) + "]"

    # split into top-level sections, collect keys per section, render, append lists
    lines = md.split("\n")
    sections, cur = [], []
    in_code = False
    for L in lines:
        if L.startswith("```"):
            in_code = not in_code
        if not in_code and L.startswith("# ") and cur:
            sections.append(cur); cur = []
        cur.append(L)
    sections.append(cur)

    out = []
    for sec in sections:
        text = "\n".join(sec)
        keys = sorted({k for m in CITE.finditer(text) for _, k, _ in _parts(m.group(1))}, key=lambda k: order[k])
        rendered = CITE.sub(render, text)
        is_biblio = "{{BIBLIOGRAPHY}}" in text
        if keys and not is_biblio:
            rendered = rendered.rstrip() + f"\n\n## {chapter_heading}\n\n" + \
                "\n\n".join(f"[{order[k]}] {refs[k][1]}" for k in keys) + "\n"
        out.append(rendered)
    md = "\n".join(out)

    full = "\n\n".join(f"[{n}] {refs[k][1]}" for k, n in sorted(order.items(), key=lambda kv: kv[1]))
    uncited = [k for k in refs if k not in order]
    if uncited:
        full += "\n\n## Further reading (not cited in the text)\n\n" + "\n\n".join(f"* {refs[k][1]}" for k in uncited)
    return md.replace("{{BIBLIOGRAPHY}}", full), order
