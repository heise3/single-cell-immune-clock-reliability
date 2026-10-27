#!/usr/bin/env python3
"""Public, exact-version source obs + Census keys. Never infer specimen from donor."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import fsspec,h5py,json,hashlib,time
from contextlib import nullcontext
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1/cache';SEED=20261025

def readcol(o,k):
 z=o[k]
 if isinstance(z,h5py.Group):
  cat=z['categories'][:];code=z['codes'][:];cat=cat.astype(str) if cat.dtype.kind in 'OS' else cat
  if code.min()<0:return np.array([cat[i] if i>=0 else '' for i in code])
  return cat[code]
 a=z[:];return a.astype(str) if a.dtype.kind in 'OS' else a

def run(r):
 import cellxgene_census as c
 out=W/f'HIHA_selected_cells_{r.dataset_id}.parquet';qa=O/f'HIHA_metadata_QA_{r.dataset_id}.json'
 if out.exists() and qa.exists():return json.loads(qa.read_text())
 cp=W/f'HIHA_Census_obs_{r.dataset_id}.parquet'
 if not cp.exists():
  if (W/'HIHA_first_Census_obs.parquet').exists() and r.dataset_id=='0a8e3443-c3e2-4918-84db-0495657d9175':cp.write_bytes((W/'HIHA_first_Census_obs.parquet').read_bytes())
  else:
   with c.open_soma(census_version='2025-11-08') as f:
    co=c.get_obs(f,'Homo sapiens',value_filter=f"dataset_id == '{r.dataset_id}'",column_names=['soma_joinid','observation_joinid','donor_id','raw_sum','cell_type']);co.to_parquet(cp,index=False)
 co=pd.read_parquet(cp).sort_values('soma_joinid').reset_index(drop=True)
 k=pd.read_csv(P/'outputs/round2/K0_primary_with_released_age.csv');kk=k[k.dataset_id.eq(r.dataset_id)].copy();want=set(kk['specimen.specimenGuid'])
 url='https://datasets.cellxgene.cziscience.com/'+r.dataset_version_id+'.h5ad';print('BEGIN',r.dataset_id,len(co),len(want),flush=True)
 local=P/'work/round2_A1/HIHA_metadata_segments'/(r.dataset_id+'.METADATA_ONLY.sparse.h5ad')
 deadline=time.monotonic()+3600
 while not local.with_suffix('.json').exists():
  if time.monotonic()>deadline:raise TimeoutError('Metadata segment not complete; no prediction permitted')
  time.sleep(10)
 with nullcontext(str(local)) as f,h5py.File(f,'r') as h:
  o=h['obs'];sample=readcol(o,'specimen.specimenGuid');sel=np.flatnonzero(np.isin(sample,list(want)));print('SAMPLE',r.dataset_id,len(sel),flush=True)
  # Compare ALL row keys, not just a sampled positional check.
  keys=readcol(o,'observation_joinid');assert len(keys)==len(co) and len(np.unique(keys))==len(keys)
  same=np.array_equal(keys,co.observation_joinid.astype(str).to_numpy())
  if same:base=co.iloc[sel].reset_index(drop=True)
  else:
   base=pd.DataFrame({'observation_joinid':keys[sel]}).merge(co,on='observation_joinid',validate='one_to_one',how='left');assert base.soma_joinid.notna().all()
  base['source_row']=sel;base['specimen']=sample[sel]
  for key in ['barcodes','subject.subjectGuid','sample.subjectAgeAtDraw','batch_id','AIFI_L1','AIFI_L2','AIFI_L3','n_umis','n_genes','pct_counts_mito']:
   base[key]=readcol(o,key)[sel];print('COL',r.dataset_id,key,flush=True)
 base['dataset_id']=r.dataset_id
 base=base.merge(kk[['gsm','donor','specimen.specimenGuid','sample.subjectAgeAtDraw','batch id']].rename(columns={'sample.subjectAgeAtDraw':'expected_age','batch id':'expected_batch'}),left_on='specimen',right_on='specimen.specimenGuid',validate='many_to_one')
 assert base.donor.eq(base['subject.subjectGuid']).all() and base.donor.eq(base.donor_id).all()
 assert np.array_equal(base['sample.subjectAgeAtDraw'].astype(int),base.expected_age.astype(int))
 # Source batch_id includes pool suffix in some releases. Preserve actual source
 # and always take batch random effect from the frozen K0 ledger.
 assert base.batch_id.astype(str).str.extract(r'^(B[0-9]+)',expand=False).eq(base.expected_batch).all()
 assert set(base.gsm)==set(kk.gsm)
 base.to_parquet(out,index=False,compression='zstd')
 q={'dataset_id':r.dataset_id,'dataset_version_id':r.dataset_version_id,'source_url':url,'n_source_cells':len(co),'n_selected_cells':len(base),'n_samples':base.gsm.nunique(),'source_Census_all_joinids_equal_in_order':bool(same),'join_method':'exact observation_joinid; full uniqueness and coverage checks','all_frozen_samples_present':True,'donor_and_age_exact_match':True,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
 qa.write_text(json.dumps(q,indent=2));print('DONE',q,flush=True);return q

def main():
 rec=pd.read_csv(O/'HIHA_Census_dataset_records.csv');results=[]
 with ThreadPoolExecutor(4) as ex:
  fs=[ex.submit(run,r) for r in rec.itertuples()]
  for f in as_completed(fs):results.append(f.result());(O/'HIHA_metadata_manifest.json').write_text(json.dumps(results,indent=2))
 a=pd.concat([pd.read_parquet(W/f'HIHA_selected_cells_{r.dataset_id}.parquet') for r in rec.itertuples()],ignore_index=True)
 assert a.gsm.nunique()==320 and a.donor.nunique()==68 and a.soma_joinid.is_unique
 a.to_parquet(W/'HIHA_K0_cells.parquet',index=False,compression='zstd');a.groupby(['donor','gsm','AIFI_L1','AIFI_L2','AIFI_L3'],observed=True).size().rename('n_cells').reset_index().to_csv(O/'K0_cell_counts_by_AIFI.csv',index=False)
 print('ALL_METADATA_DONE',len(a),flush=True)
if __name__=='__main__':main()
