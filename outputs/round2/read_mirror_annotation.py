#!/usr/bin/env python3
"""Read only public CELLxGENE var/uns annotation via HTTP range requests."""
from pathlib import Path
import json
import fsspec,h5py,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent;PROJECT=ROOT.parents[1];WORK=PROJECT/'work/round2'
assets=json.loads((WORK/'K0_snapshot/K0_CELLxGENE_assets.json').read_text())['datasets']
def read_col(obj):
 if isinstance(obj,h5py.Group):
  cats=obj['categories'][:];codes=obj['codes'][:];return np.array([cats[i] if i>=0 else '' for i in codes])
 return obj[:]
records=[]
# Each dataset's var is checked independently; no assumption that all mirrors match.
for d in assets:
 u=next(a['url'] for a in d['assets'] if a['filetype']=='H5AD')
 with fsspec.open(u,'rb',block_size=2**20,cache_type='blockcache',client_kwargs={'trust_env':True}) as stream:
  with h5py.File(stream,'r') as f:
   var=f['var'];data={k:read_col(var[k]) for k in var.keys()}
   for k,v in data.items():
    if v.dtype.kind in ['O','S']:data[k]=v.astype(str)
   pd.DataFrame(data).to_csv(ROOT/f'CELLxGENE_var_{d["dataset_id"]}.csv',index=False)
   uns={}
   if 'uns' in f:
    for k in f['uns']:
     if any(t in k.lower() for t in ['schema','citation','title','organism','version','reference']):
      o=f['uns'][k]
      if isinstance(o,h5py.Dataset):uns[k]=str(o[()])[:2000]
   records.append({'dataset_id':d['dataset_id'],'dataset_version_id':d['dataset_version_id'],'url':u,'var_columns':list(data),'n_features':len(next(iter(data.values()))),'uns':uns})
 print(records[-1],flush=True)
(ROOT/'CELLxGENE_annotation_inspection.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
