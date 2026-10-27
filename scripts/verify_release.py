from pathlib import Path
import ast, hashlib, json
P=Path(__file__).resolve().parents[1]
rows=(P/'SHA256SUMS.txt').read_text().splitlines()
for line in rows:
    digest,name=line.split('  ',1)
    assert hashlib.sha256((P/name).read_bytes()).hexdigest()==digest,name
files=list(P.rglob('*.py'))
for p in files:
    if '.git' not in p.parts:ast.parse(p.read_text(),filename=str(p))
print(json.dumps({'hashes_verified':len(rows),'python_sources_parsed':len(files),'empirical_reanalysis':False}))
