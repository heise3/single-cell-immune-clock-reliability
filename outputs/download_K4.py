#!/usr/bin/env python3
"""Download only the frozen 20-donor/40-visit public GEO K4 files."""
import json,hashlib,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parent
DATA=ROOT.parent/'work/K4_data';DATA.mkdir(exist_ok=True)
plan=json.loads((ROOT/'K4_download_plan.json').read_text())
def fetch(item):
 p=DATA/item['filename'];tmp=p.with_suffix('.partial')
 if not p.exists():
  r=subprocess.run(['curl','-fsSL','--max-time','150','--retry','1','-o',str(tmp),item['url']],capture_output=True,text=True)
  if r.returncode:return {**item,'error':r.stderr.strip()}
  tmp.rename(p)
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return {**item,'bytes':p.stat().st_size,'sha256':h.hexdigest(),'path':str(p.relative_to(ROOT.parent))}
results=[]
with ThreadPoolExecutor(max_workers=4) as pool:
 futures=[pool.submit(fetch,item) for item in plan['files']]
 for future in as_completed(futures):
  r=future.result();results.append(r);print(len(results),r['gsm'],r.get('bytes',r.get('error')),flush=True)
(ROOT/'K4_download_manifest.json').write_text(json.dumps({'seed':20261025,'files':sorted(results,key=lambda r:r['gsm'])},indent=2))
if any('error' in r for r in results):raise SystemExit('Download failures: no outcome analysis authorized on incomplete data.')
