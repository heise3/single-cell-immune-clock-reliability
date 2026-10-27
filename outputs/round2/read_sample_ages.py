#!/usr/bin/env python3
"""Read released subjectAgeAtDraw and join by sample ID; do not invent exact age."""
from pathlib import Path
import json
import fsspec,h5py,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent;PROJECT=ROOT.parents[1];WORK=PROJECT/'work/round2'
assets=json.loads((WORK/'K0_snapshot/K0_CELLxGENE_assets.json').read_text())['datasets']
def col(o):
 if isinstance(o,h5py.Group):
  a=o['categories'][:];codes=o['codes'][:];return np.array([a[i] if i>=0 else '' for i in codes])
 return o[:]
all_rows=[]
for d in assets:
 cache=WORK/'cache'/f'sample_ages_{d["dataset_id"]}.csv'
 if cache.exists():all_rows.append(pd.read_csv(cache));continue
 u=next(a['url'] for a in d['assets'] if a['filetype']=='H5AD')
 with fsspec.open(u,'rb',block_size=2**20,cache_type='blockcache',client_kwargs={'trust_env':True}) as stream:
  with h5py.File(stream,'r') as f:
   keys=['sample.sampleKitGuid','specimen.specimenGuid','subject.subjectGuid','sample.subjectAgeAtDraw','subject.ageAtFirstDraw','subject.birthYear','sample.drawYear']
   values={k:col(f['obs'][k]) for k in keys}
   for k in values:
    if values[k].dtype.kind in ['S','O']:values[k]=values[k].astype(str)
   data=pd.DataFrame(values)
   assert data.groupby('sample.sampleKitGuid')['sample.subjectAgeAtDraw'].nunique().max()==1
   data=data.drop_duplicates();data['dataset_id']=d['dataset_id'];data.to_csv(cache,index=False);all_rows.append(data)
 print('ages',d['dataset_id'],len(data),flush=True)
out=pd.concat(all_rows,ignore_index=True);out.to_csv(ROOT/'released_sample_age_metadata.csv',index=False)
main=pd.read_csv(ROOT/'K0_primary_68_donors_320_visits.csv')
joined=main.merge(out,on=['specimen.specimenGuid','dataset_id'],validate='one_to_one')
assert len(joined)==320 and joined['sample.subjectAgeAtDraw'].notna().all()
joined.to_csv(ROOT/'K0_primary_with_released_age.csv',index=False)
