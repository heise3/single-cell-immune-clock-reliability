#!/usr/bin/env python3
"""Acquire only frozen K0 68/320 public matrices, preserving original caches."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import requests,pandas as pd,numpy as np,h5py,json,hashlib,os,time
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';D=P/'work/round2_A1/main_GEO';D.mkdir(exist_ok=True)
def main():
 g=pd.read_csv(P/'outputs/round2/K0_primary_with_released_age.csv');m=pd.read_csv(P/'outputs/sample_metadata.csv')
 plan=g.merge(m[['gsm','count_url']],on='gsm',validate='one_to_one');assert len(plan)==320 and plan.donor.nunique()==68
 plan['filename']=plan.count_url.str.split('/').str[-1];plan['expected_pbmc_sample_id']=plan['specimen.specimenGuid'];plan.to_csv(O/'main_GEO_download_plan.csv',index=False)
 reference=pd.read_csv(P/'outputs/round2/GEO_gene_reference.csv');expected=reference.geo_feature_id.to_numpy();results=[]
 def get(row):
  path=D/row.filename;old=P/'work/K4_data'/row.filename;rec={'gsm':row.gsm,'donor':row.donor,'filename':row.filename,'url':row.count_url,'path':str(path.relative_to(P))}
  try:
   if not path.exists() and old.exists():os.link(old,path);rec['route']='reuse_immutable_pilot_cache'
   if not path.exists():
    rec['route']='anonymous_GEO_HTTPS';part=path.with_suffix(path.suffix+'.part')
    for attempt in range(3):
     try:
      with requests.get(row.count_url,stream=True,timeout=(15,90)) as r:
       r.raise_for_status()
       with part.open('wb') as f:
        for b in r.iter_content(2**20):f.write(b)
      with h5py.File(part,'r') as f:assert 'matrix' in f
      part.rename(path);break
     except Exception:
      if attempt==2:raise
   with h5py.File(path,'r') as f:
    z=f['matrix'];ids=z['features/id'][:].astype(str);rec['n_genes'],rec['n_cells']=map(int,z['shape'][:]);rec['features_match_frozen_reference']=bool(np.array_equal(ids,expected))
    sample=np.unique(z['observations/pbmc_sample_id'][:].astype(str)).tolist();assert sample==[row.expected_pbmc_sample_id]
    rec['pbmc_sample_id']=sample[0];assert sample[0] in row.filename
   rec['bytes']=path.stat().st_size;rec['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();rec['status']='DOWNLOADED_VERIFIED'
  except Exception as e:rec['status']='ACQUISITION_FAILED';rec['error']=repr(e)
  return rec
 with ThreadPoolExecutor(6) as ex:
  fut=[ex.submit(get,r) for r in plan.itertuples()]
  for f in as_completed(fut):
   rec=f.result();results.append(rec);print(len(results),rec['gsm'],rec['status'],rec.get('bytes'),flush=True)
   (O/'main_GEO_download_manifest.json').write_text(json.dumps({'n_expected':320,'files':results},indent=2))
 assert all(r['status']=='DOWNLOADED_VERIFIED' for r in results),'Acquisition failure, not a scientific KILL.'
 assert all(r['features_match_frozen_reference'] for r in results),'Feature changes require explicit sample-specific mapping.'
 print('ALL_DONE_320',flush=True)
if __name__=='__main__':main()
