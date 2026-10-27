#!/usr/bin/env python3
from pathlib import Path
import h5py,fsspec,aiohttp,pandas as pd,numpy as np,json,re,time,sys
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';W=P/'work/round3'
def col(o):
 if isinstance(o,h5py.Group):
  codes=o['codes'][:];cats=o['categories'][:];return np.array([cats[i] if i>=0 else '' for i in codes]).astype(str)
 a=o[:];return a.astype(str) if a.dtype.kind in 'OS' else a
def main():
 r=pd.read_csv(O/'AIDA_Census_selected_v2.csv').iloc[0];url=f'https://datasets.cellxgene.cziscience.com/{r.dataset_version_id}.h5ad';fs=fsspec.filesystem('https',client_kwargs={'trust_env':True,'timeout':aiohttp.ClientTimeout(total=600)});report={'url':url,'dataset_id':r.dataset_id,'dataset_version_id':r.dataset_version_id,'purpose':'metadata only; never infer raw count contents from this read'}
 sys.path.insert(0,str(P/'outputs/round2_A1'));import cache_HIHA_annotation_segments as cache
 cache.D=W/'AIDA_metadata_segments';cache.D.mkdir(exist_ok=True);receipt=cache.run(r);report['metadata_segment_receipt']=receipt;local=cache.D/(r.dataset_id+'.METADATA_ONLY.sparse.h5ad')
 from metadata_range_reader import MetadataRangeReader
 reader=MetadataRangeReader(local,url,fs,receipt,local.with_suffix('.json'))
 with reader,h5py.File(reader,'r') as h:
  keys=list(h['obs']);report['obs_columns']=keys;print('SCHEMA',keys,flush=True);d={}
  for k in keys:
   if k in ['development_stage','Country','donor_id','self_reported_ethnicity','sex']:d[k]=col(h['obs'][k]);print('COLUMN',k,len(np.unique(d[k])),np.unique(d[k])[:10],flush=True)
  pd.DataFrame(d).to_parquet(W/'cache/AIDA_source_obs.parquet',index=False)
  v=pd.DataFrame({k:col(z) for k,z in h['var'].items()});v.to_csv(O/'AIDA_source_var.csv',index=False);report['n_source_features']=len(v);report['n_source_cells']=len(next(iter(d.values())))
 (O/'AIDA_source_metadata_preflight.json').write_text(json.dumps(report,indent=2));print('DONE',flush=True)
if __name__=='__main__':main()
