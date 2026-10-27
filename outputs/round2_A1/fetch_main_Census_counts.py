#!/usr/bin/env python3
"""Frozen K0 exact keys, raw counts, full-library normalization, resumable chunks."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';D=W/'main_Census_counts';D.mkdir(exist_ok=True)
L2={'Memory CD4 T cell':'CD4T','Naive CD4 T cell':'CD4T','Treg':'CD4T','Memory CD8 T cell':'CD8T','Naive CD8 T cell':'CD8T','Naive B cell':'B','Memory B cell':'B','Transitional B cell':'B','Effector B cell':'B','CD14 monocyte':'MONO','CD16 monocyte':'MONO','Intermediate monocyte':'MONO','CD56bright NK cell':'NK','CD56dim NK cell':'NK','Proliferating NK cell':'NK'}
def main():
 import cellxgene_census as c
 manifest_path=O/'main_Census_count_manifest.json'
 if manifest_path.exists():
  cached=json.loads(manifest_path.read_text())
  if cached.get('complete') and cached.get('exact_union_storage_reuse'):
   assert sum(r['n_cells'] for r in cached['chunks'])==cached['plan']['n_cells']
   for r in cached['chunks']:
    source=P/r['path'];assert source.stat().st_size==r['bytes'];assert hashlib.sha256(source.read_bytes()).hexdigest()==r['sha256']
   print('VERIFIED_FROZEN_POPULATION_UNION_CACHE',cached['plan'],flush=True);return
 obs=pd.read_parquet(W/'cache/HIHA_K0_cells.parquet',columns=['soma_joinid','gsm','donor','raw_sum','AIFI_L2','AIFI_L3','dataset_id']);obs['clock_type']=obs.AIFI_L2.map(L2).fillna('')
 ram=pd.read_csv(O/'RAM_mapped_gene_panels.csv',keep_default_na=False);rtypes=set(ram.RAM_celltype)
 # Exclude ambiguous T lineages, plasma cells and non-clock lineages; no quality
 # thresholds re-estimated after outcomes. Keep original released QC cells.
 obs=obs[obs.clock_type.ne('')|obs.AIFI_L3.isin(rtypes)].sort_values('soma_joinid').reset_index(drop=True)
 assert obs.gsm.nunique()==320 and obs.donor.nunique()==68;obs.to_parquet(W/'cache/HIHA_K0_clock_cells.parquet',index=False)
 co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');m=pd.read_csv(O/'gene_mapping_HIHA_Census.csv',keep_default_na=False);targets=set(co.loc[co.gene.ne('(Intercept)')&co.coefficient.ne(0),'gene']);genes={e for s in m.loc[m.model_ensembl.isin(targets),'mapped_feature_ids'] for e in s.split('|') if e};genes|={e for s in ram.mapped_feature_ids for e in s.split('|') if e}
 v=pd.read_csv(O/'HIHA_Census_var_0a8e3443-c3e2-4918-84db-0495657d9175.csv');v=v[v['_index'].isin(genes)].sort_values('soma_joinid');assert len(v)==len(genes)
 plan={'soma_buffer_bytes':67108864,'TileDB_total_budget_bytes':8589934592,'seed':20261025,'Census_release':'2025-11-08','n_cells':len(obs),'n_features':len(v),'n_donors':68,'n_visits':320,'label_crosswalk':L2,'source_QC':'all final source-released cells in eligible lineages; no additional outcome-based QC','library_denominator':'Census raw_sum full released feature universe','exact_identity_precheck':'all source observation_joinids matched one-to-one to Census'};(O/'main_Census_execution_plan.json').write_text(json.dumps(plan,indent=2));print('PLAN',plan,flush=True)
 rows=[]
 with c.open_soma(census_version='2025-11-08',tiledb_config={'soma.init_buffer_bytes':67108864,'py.init_buffer_bytes':67108864,'sm.mem.total_budget':8589934592,'sm.compute_concurrency_level':4,'sm.io_concurrency_level':4}) as f:
  for i,lo in enumerate(range(0,len(obs),100000)):
   rr=obs.iloc[lo:lo+100000];p=D/f'HIHA_K0_{i:03d}.h5ad';print('QUERY',i,len(rr),flush=True)
   if not p.exists():
    a=c.get_anndata(f,'Homo sapiens',obs_coords=rr.soma_joinid.to_numpy(dtype=np.int64),var_coords=v.soma_joinid.to_numpy(dtype=np.int64),obs_column_names=['soma_joinid','raw_sum'],var_column_names=['soma_joinid','feature_id','feature_name'])
    assert a.n_obs==len(rr) and a.n_vars==len(v);a.obs=a.obs.reset_index(drop=True).merge(rr.drop(columns=['raw_sum']),on='soma_joinid',validate='one_to_one');a.obs.index=a.obs.soma_joinid.astype(str);a.obs.index.name=None;assert np.all(a.X.data>=0) and np.allclose(a.X.data,np.rint(a.X.data));assert np.all(np.asarray(a.X.sum(axis=1)).ravel()<=a.obs.raw_sum.to_numpy()+1e-5)
    a.write_h5ad(p.with_suffix('.part.h5ad'),compression='gzip');p.with_suffix('.part.h5ad').rename(p);del a
   rows.append({'path':str(p.relative_to(P)),'n_cells':len(rr),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});(O/'main_Census_count_manifest.json').write_text(json.dumps({'plan':plan,'chunks':rows},indent=2));print('DONE',i,flush=True)
 print('ALL_RAW_COUNTS_DONE',flush=True)
if __name__=='__main__':main()
