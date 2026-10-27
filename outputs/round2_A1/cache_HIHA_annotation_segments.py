#!/usr/bin/env python3
"""Sparse local metadata cache: header plus complete source obs-to-EOF segment.
ZERO HOLE IS NOT EXPRESSION. Never read X from these files. Exact URL pinned.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import fsspec,aiohttp,h5py,json,hashlib,time
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';D=P/'work/round2_A1/HIHA_metadata_segments';D.mkdir(exist_ok=True)
def run(r):
 p=D/(r.dataset_id+'.METADATA_ONLY.sparse.h5ad');j=p.with_suffix('.json')
 if p.exists() and j.exists():return json.loads(j.read_text())
 url='https://datasets.cellxgene.cziscience.com/'+r.dataset_version_id+'.h5ad';fs=fsspec.filesystem('https',client_kwargs={'trust_env':True,'timeout':aiohttp.ClientTimeout(total=600)})
 info=fs.info(url);size=info['size'];head=fs.cat_file(url,start=0,end=1048576)
 with p.open('wb') as f:f.truncate(size);f.write(head)
 with fs.open(url,'rb',block_size=1048576,cache_type='blockcache') as remote,h5py.File(remote,'r') as h:
  root=h['/'];start=int(root.id.links.get_info(b'obs').u)
 assert 1048576<start<size
 print('TAIL_BEGIN',r.dataset_id,start,size,size-start,flush=True);t=time.monotonic()
 # Contiguous 32MB ranges; public client handles HTTP, not raw network access.
 hashes=[]
 with p.open('r+b') as f:
  for lo in range(start,size,32*1024**2):
   hi=min(lo+32*1024**2,size);b=fs.cat_file(url,start=lo,end=hi);assert len(b)==hi-lo,(len(b),hi-lo);f.seek(lo);f.write(b);hashes.append({'start':lo,'end_exclusive':hi,'sha256':hashlib.sha256(b).hexdigest()});print('TAIL',r.dataset_id,hi-start,'/',size-start,flush=True)
 with h5py.File(p,'r') as h:assert 'observation_joinid' in h['obs']
 z={'dataset_id':r.dataset_id,'dataset_version_id':r.dataset_version_id,'url':url,'logical_file_size':size,'copied_segments':[[0,1048576],[start,size]],'zero_hole':[1048576,start],'NOT_A_FULL_EXPRESSION_FILE':True,'header_sha256':hashlib.sha256(head).hexdigest(),'tail_chunks':hashes,'seconds':time.monotonic()-t};j.write_text(json.dumps(z,indent=2));return z
if __name__=='__main__':
 rows=pd.read_csv(O/'HIHA_Census_dataset_records.csv');out=[]
 with ThreadPoolExecutor(4) as ex:
  for z in as_completed([ex.submit(run,r) for r in rows.itertuples()]):out.append(z.result());(O/'HIHA_annotation_segment_manifest.json').write_text(json.dumps(out,indent=2))
