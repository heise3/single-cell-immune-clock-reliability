#!/usr/bin/env python3
"""Refresh only document corrections and the archive; do not rerun scientific QA."""
from pathlib import Path
import sys,json,zipfile,hashlib,datetime
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 qa=json.loads((O/'final_QA.json').read_text());n=qa['count_cache_files_verified_unchanged'];f=O/'README_reproduce.md';s=f.read_text().replace('142个count缓存',f'{n}个count缓存');f.write_text(s);f=O/'package_round5.py';f.write_text(f.read_text().replace('142个count缓存',f'{n}个count缓存'))
 note={'scope':'Documentation-only correction:142count files typo corrected to verified143; science inputs/results unchanged','source_note':'Initial unused PALMO OA-ID request404 was superseded by correct primary PMC10041512 XML;not interpreted as data absence','scientific_QA':'prior final_QA retained,not recomputed','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(O/'documentation_finalize_receipt.json').write_text(json.dumps(note,indent=2))
 rows=[]
 for p in sorted(O.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS.csv':rows.append({'file':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':sha(p)})
 pd.DataFrame(rows).to_csv(O/'SHA256SUMS.csv',index=False);target=P/'outputs/round5_checkpoint.zip';tmp=target.with_suffix('.part.zip')
 with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for r in rows:z.write(O/r['file'],arcname='round5/'+r['file'])
  z.write(O/'SHA256SUMS.csv',arcname='round5/SHA256SUMS.csv')
 tmp.replace(target)
 with zipfile.ZipFile(target) as z:
  assert z.testzip() is None
  for r in rows:assert hashlib.sha256(z.read('round5/'+r['file'])).hexdigest()==r['sha256']
 receipt={'archive':str(target.relative_to(P)),'bytes':target.stat().st_size,'sha256':sha(target),'manifest_files':len(rows),'CRC_verified':True,'all_archive_member_SHA256_verified':True,'seed':20261025,'status':'COMPLETE','preregistration_P2_label':'待声明','historical_archives_unchanged':True,'recorded_UTC':note['recorded_UTC']};(P/'outputs/round5_checkpoint_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2));print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
