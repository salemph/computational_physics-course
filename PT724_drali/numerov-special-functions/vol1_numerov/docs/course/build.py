import re, subprocess, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / 'course'))
from refs_numerov import REFS
from cite import process
s = pathlib.Path('course.md').read_text()
def rep(m):
    lang, path = m.group(1), m.group(2)
    code = pathlib.Path(path).read_text().rstrip()
    name = pathlib.Path(path).name
    return f"File: `{name}`\n\n```{lang}\n{code}\n```"
s, ORDER = process(s, REFS, further=True)
print(len(ORDER), 'references cited')
s = re.sub(r"\{\{CODE:(\w+):([^}]+)\}\}", rep, s)
pathlib.Path('course_full.md').write_text(s)
subprocess.run(["pandoc", "course_full.md", "-o", "raw.docx", "--toc", "--toc-depth=2",
                "--highlight-style=tango", "--resource-path=.", "-f", "markdown+tex_math_dollars"], check=True)
print("raw.docx built")
