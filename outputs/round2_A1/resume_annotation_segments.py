#!/usr/bin/env python3
"""Resume only explicitly logged committed ranges; parallel immutable HTTP reads."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import re,json,hashlib,time
import fsspec,aiohttp,pandas as pd,h5py
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';D=W/'HIHA_metadata_segments';log=(W/'HIHA_annotation_segments_retry.log').read_text();start={};done={}
for line in log.splitlines():
 a=line.split()
 if a and a[0]=='TAIL_BEGIN':start[a[1]]=(int(a[2]),int(a[3]))
 if a and a[0]=='TAIL':done[a[1]]=max(done.get(a[1],0),int(a[2]))
records=pd.read_csv(O/'HIHA_Census_dataset_records.csv');state={};jobs=[];fs=fsspec.filesystem('https',client_kwargs={'trust_env':True,'timeout':aiohttp.ClientTimeout(total=240)})
for r in records.itertuples():
 p=D/(r.dataset_id+'.METADATA_ONLY.sparse.h5ad');j=p.with_suffix('.json')
 if j.exists():state[r.dataset_id]=json.loads(j.read_text());continue
 lo,size=start[r.dataset_id];assert p.stat().st_size==size;url='https://datasets.cellxgene.cziscience.com/'+r.dataset_version_id+'.h5ad';chunks=[]
 with p.open('rb') as f:
  head=f.read(1048576)
  for off in range(lo,size,32*1024**2):
   hi=min(off+32*1024**2,size)
   if hi-lo<=done.get(r.dataset_id,0):f.seek(off);b=f.read(hi-off);chunks.append({'start':off,'end_exclusive':hi,'sha256':hashlib.sha256(b).hexdigest(),'route':'verified formerly logged complete HTTP range'})
   else:jobs.append((r.dataset_id,str(p),url,off,hi))
 state[r.dataset_id]={'dataset_id':r.dataset_id,'dataset_version_id':r.dataset_version_id,'url':url,'logical_file_size':size,'copied_segments':[[0,1048576],[lo,size]],'zero_hole':[1048576,lo],'NOT_A_FULL_EXPRESSION_FILE':True,'header_sha256':hashlib.sha256(head).hexdigest(),'tail_chunks':chunks,'expected_tail_chunks':len(range(lo,size,32*1024**2))}

def fetch(job):
 did,path,url,lo,hi=job
 for retry in range(3):
  try:
   b=fs.cat_file(url,start=lo,end=hi);assert len(b)==hi-lo
   with open(path,'r+b') as f:f.seek(lo);f.write(b);f.flush()
   return did,{'start':lo,'end_exclusive':hi,'sha256':hashlib.sha256(b).hexdigest(),'route':'public HTTP exact range'}
  except Exception:
   if retry==2:raise
 print('UNREACHABLE')
print('PENDING_RANGES',len(jobs),flush=True)
with ThreadPoolExecutor(8) as ex:
 for future in as_completed([ex.submit(fetch,job) for job in jobs]):
  did,chunk=future.result();z=state[did];z['tail_chunks'].append(chunk);print('RANGE_DONE',did,chunk['start'],len(z['tail_chunks']),'/',z['expected_tail_chunks'],flush=True)
  p=D/(did+'.METADATA_ONLY.sparse.h5ad');p.with_suffix('.progress.json').write_text(json.dumps(z,indent=2))
  if len(z['tail_chunks'])==z['expected_tail_chunks']:
   z['tail_chunks']=sorted(z['tail_chunks'],key=lambda x:x['start'])
   with h5py.File(p,'r') as h:assert 'observation_joinid' in h['obs']
   p.with_suffix('.json').write_text(json.dumps(z,indent=2));print('DATASET_COMPLETE',did,flush=True)
(O/'HIHA_annotation_segment_manifest.json').write_text(json.dumps(list(state.values()),indent=2));print('ALL_SEGMENTS_COMPLETE',flush=True)
