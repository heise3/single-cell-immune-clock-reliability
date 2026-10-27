#!/usr/bin/env python3
"""Public range-only preflight. Never treat metadata availability as matrix QA."""
from pathlib import Path
import argparse,json,re
import fsspec,h5py,numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1/cache'
def col(o):
 if isinstance(o,h5py.Group):
  if 'codes' in o and 'categories' in o:
   c=o['categories'][:];i=o['codes'][:];return np.array([c[k] if k>=0 else '' for k in i]).astype(str)
  raise ValueError(('unsupported H5 column',o.name,list(o.keys())))
 a=o[:];return a.astype(str) if a.dtype.kind in 'OS' else a
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--metadata');ap.add_argument('--url');ap.add_argument('--stem',required=True);a=ap.parse_args()
 m=json.loads(Path(a.metadata).read_text()) if a.metadata else {}
 u=a.url or next(z['url'] for z in m['assets'] if z['filetype']=='H5AD')
 d={'url':u,'title':m.get('title'),'dataset_version_id':m.get('dataset_version_id'),'cell_count_from_catalog':m.get('cell_count')}
 with fsspec.open(u,'rb',block_size=2**20,cache_type='blockcache',client_kwargs={'trust_env':True}) as s:
  with h5py.File(s,'r') as f:
   d['obs_columns']=list(f['obs']);d['root_keys']=list(f);print('SCHEMA',d,flush=True)
   for key in ['var','raw/var']:
    if key not in f:continue
    v=pd.DataFrame({k:col(z) for k,z in f[key].items()});v.to_csv(O/(a.stem+'_'+key.replace('/','_')+'.csv'),index=False)
    d[key]={'n_features':len(v),'columns':list(v)};print(key,d[key],flush=True)
   for k in ['X','raw/X']:
    if k in f:
     z=f[k];d[k]={'shape':list(z.attrs.get('shape',z.shape if isinstance(z,h5py.Dataset) else [])),'encoding':str(z.attrs.get('encoding-type','')),'children':list(z) if isinstance(z,h5py.Group) else None}
   selected=[k for k in f['obs'] if re.search(r'(^age$|age_|_age|donor|individual|subject|cell.*type|cell.*label|annotation|^sex$|^disease$|^tissue$|sample)',k,re.I)]
   obs={}
   for k in selected:
    obs[k]=col(f['obs'][k]);print('OBS',k,'unique',len(np.unique(obs[k])),'sample',np.unique(obs[k])[:8].tolist(),flush=True)
   pd.DataFrame(obs).to_csv(W/(a.stem+'_selected_obs.csv.gz'),index=False,compression='gzip')
 d={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in d.items()}
 (O/(a.stem+'_preflight.json')).write_text(json.dumps(d,indent=2,default=int));print('DONE',a.stem,flush=True)
if __name__=='__main__':main()
