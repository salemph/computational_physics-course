import re, subprocess, pathlib, sys
from refs import REFS
from cite import process
s = pathlib.Path('course.md').read_text()
def rep(m):
    lang, path = m.group(1), m.group(2)
    code = pathlib.Path(path).read_text().rstrip()
    name = pathlib.Path(path).name
    return f"File: `{name}`\n\n```{lang}\n{code}\n```"
s, ORDER = process(s, REFS)
s = re.sub(r"\{\{CODE:(\w+):([^}]+)\}\}", rep, s)
print(len(ORDER), 'references cited')
# proportional column widths: rewrite each table separator line from the column contents
lines = s.split('\n'); out = []
i = 0
while i < len(lines):
    L = lines[i]
    if i + 1 < len(lines) and L.startswith('|') and re.match(r'^\|(\s*:?-+:?\s*\|)+\s*$', lines[i + 1]):
        j = i + 2
        while j < len(lines) and lines[j].startswith('|'): j += 1
        rows = [lines[i]] + lines[i + 2:j]
        splitp = lambda r: [c.strip() for c in re.split(r'(?<!\\)\|', r.strip()[1:-1])]
        cells = [splitp(r) for r in rows]
        ncol = len(cells[0])
        w = [max(len(re.sub(r'[\$\\{}*]', '', r[k])) if k < len(r) else 1 for r in cells) for k in range(ncol)]
        seps = [c.strip() for c in lines[i + 1].strip().strip('|').split('|')]
        new = []
        for k in range(ncol):
            right = seps[k].endswith(':') if k < len(seps) else False
            d = '-' * max(7, w[k])
            new.append(d[:-1] + ':' if right else d)
        out += [lines[i], '|' + '|'.join(new) + '|'] + lines[i + 2:j]
        i = j
    else:
        out.append(L); i += 1
s = '\n'.join(out)
pathlib.Path('course_full.md').write_text(s)
subprocess.run(["pandoc", "course_full.md", "-o", "raw.docx", "--toc", "--toc-depth=2",
                "--highlight-style=tango", "--resource-path=.", "-f", "markdown+tex_math_dollars"], check=True)
print("raw.docx built")
