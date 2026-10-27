#!/usr/bin/env python3
"""Download only the registered source after all K gates; fixed SHA objects."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,time,requests,hashlib,datetime,tarfile,shutil
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';W=P/'work/round6';D=W/'source';E=W/'extracted'
def one(item):
 ct,r=item;name=r['name'];target=D/name;url='https://humancellatlas.s3.amazonaws.com/'+r['azul_mirror_uri'].split('s3://humancellatlas/',1)[1];z={'celltype':ct,'url':url,'expected_sha256':r['sha256'],'expected_bytes':r['size'],'start_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'anonymous':True};h=hashlib.sha256();n=0;start=time.monotonic()
 try:
  if not target.exists():
   with requests.get(url,stream=True,timeout=(45,90)) as resp:
    z['http_status']=resp.status_code;resp.raise_for_status()
    with target.with_suffix('.gz.part').open('wb') as f:
     for b in resp.iter_content(4*1024**2):
      f.write(b);h.update(b);n+=len(b)
      if n//134217728!=(n-len(b))//134217728:print('MATRIX_BYTES',ct,n,'/',r['size'],flush=True)
      if time.monotonic()-start>1800:raise TimeoutError('1800sec individualobjectcap')
   assert n==r['size'] and h.hexdigest()==r['sha256'],(name,n,h.hexdigest());target.with_suffix('.gz.part').replace(target)
  else:
   with target.open('rb') as f:
    for b in iter(lambda:f.read(8388608),b''):h.update(b);n+=len(b)
   assert n==r['size'] and h.hexdigest()==r['sha256']
  z.update(bytes=n,sha256=h.hexdigest(),path=str(target.relative_to(P)));members=[]
  with tarfile.open(target) as t:
   for m in t:
    members.append({'name':m.name,'size':m.size})
    if m.isfile() and m.name.endswith('_rna.rds'):
     f=E/(ct+'_rna.rds');hs=hashlib.sha256()
     with t.extractfile(m) as s,f.open('wb') as out:
      for b in iter(lambda:s.read(8388608),b''):out.write(b);hs.update(b)
     z.update(RNA_RDS_path=str(f.relative_to(P)),RNA_RDS_sha256=hs.hexdigest(),RNA_RDS_bytes=f.stat().st_size)
  assert 'RNA_RDS_path'in z;z.update(members=members,status='COMPLETE',seconds=time.monotonic()-start)
 except Exception as e:z.update(status='FAILED',error=repr(e),seconds=time.monotonic()-start)
 return z
def main():
 assert all(json.loads((O/f'K_{g}_status.json').read_text())['status']=='PASS' for g in ['a','b','c']);refs=json.loads((O/'HCA_declared_file_records.json').read_text());jobs=[(ct,refs[name]) for ct,name in [('CD4T','cd4_t_cells.tar.gz'),('CD8T','conventional_cd8_t_cells.tar.gz'),('MONO','myeloid_cells.tar.gz'),('NK','nk_cells.tar.gz')]];rows=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(one,x) for x in jobs]):
   z=f.result();rows.append(z);(O/'primary_matrix_download_manifest.json').write_text(json.dumps(rows,indent=2));print('MATRIX_DONE',z['celltype'],z['status'],z.get('seconds'),z.get('error',''),flush=True)
 assert all(x['status']=='COMPLETE' for x in rows),rows
if __name__=='__main__':main()
