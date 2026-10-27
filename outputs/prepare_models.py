#!/usr/bin/env python3
"""Fetch the frozen public pretrained-model releases used by this audit."""
from pathlib import Path
import subprocess,zipfile,hashlib,json
ROOT=Path(__file__).resolve().parent;WORK=ROOT.parent/'work';WORK.mkdir(exist_ok=True)
items=[
 ('https://codeload.github.com/CiiM-Bioinformatics-group/scImmuAging/zip/112f8fab10c72aec6e7a7d83762a40117a1e6e24','scImmuAging.zip'),
 ('https://codeload.github.com/cliffren/PHARE/zip/0ab277b38ddac534d6d91e4580134cf6c78822ca','PHARE.zip'),
 ('https://raw.githubusercontent.com/shanzha9/siAgeModel/3e8e17567afa296213b414083eaca59b30455b81/model/all_selectedModel.rds','siAgeModel.rds')]
records=[]
for url,name in items:
 p=WORK/name
 if not p.exists():subprocess.run(['curl','-fsSL','--max-time','60','-o',str(p),url],check=True)
 if p.suffix=='.zip':
  with zipfile.ZipFile(p) as z:z.extractall(WORK/'models')
 records.append({'url':url,'path':str(p.relative_to(ROOT.parent)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(ROOT/'model_manifest.json').write_text(json.dumps(records,indent=2))
