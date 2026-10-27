#!/usr/bin/env python3
"""Fixed age-stratified 200-donor technical-mask input, public Census raw counts.
Selection never uses predictions. The original 981-donor age ledger is retained.
Queries contiguous ranges of selected soma_joinids, to avoid rereading tiles.
"""
from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';D=W/'OneK1K_counts';D.mkdir(exist_ok=True)
SEED=20261025;RELEASE='2025-11-08'
CROSSWALK={
 'central memory CD4-positive, alpha-beta T cell':'CD4T','naive thymus-derived CD4-positive, alpha-beta T cell':'CD4T',
 'effector memory CD4-positive, alpha-beta T cell':'CD4T','regulatory T cell':'CD4T','CD4-positive, alpha-beta cytotoxic T cell':'CD4T','CD4-positive, alpha-beta T cell':'CD4T',
 'effector memory CD8-positive, alpha-beta T cell':'CD8T','naive thymus-derived CD8-positive, alpha-beta T cell':'CD8T',
 'central memory CD8-positive, alpha-beta T cell':'CD8T','CD8-positive, alpha-beta T cell':'CD8T',
 'naive B cell':'B','memory B cell':'B','transitional stage B cell':'B',
 'CD14-positive monocyte':'MONO','CD14-low, CD16-positive monocyte':'MONO',
 'natural killer cell':'NK','CD16-negative, CD56-bright natural killer cell, human':'NK'}
def main():
 import cellxgene_census as c
 import anndata as ad
 obs=pd.read_parquet(W/'cache/OneK1K_census_obs.parquet');obs['broad_type']=obs.cell_type.astype(str).map(CROSSWALK)
 ages=pd.read_csv(W/'cache/OneK1K_selected_obs.csv.gz',usecols=['donor_id','age']).drop_duplicates()
 assert ages.groupby('donor_id').size().max()==1
 ages=ages.sort_values(['age','donor_id']).reset_index(drop=True);ages['age_stratum']=pd.qcut(ages.age,5,labels=False,duplicates='drop')
 rng=np.random.default_rng(SEED);chosen=[]
 for _,g in ages.groupby('age_stratum'):chosen.extend(rng.choice(g.donor_id.astype(str).to_numpy(),40,replace=False).tolist())
 # Fix both observed age endpoints in their original age strata.
 for ix in [0,len(ages)-1]:
  row=ages.iloc[ix];donor=str(row.donor_id)
  if donor not in chosen:
   candidates=ages[ages.age_stratum.eq(row.age_stratum)&ages.donor_id.isin(chosen)].sort_values('donor_id')
   chosen.remove(str(candidates.iloc[-1].donor_id));chosen.append(donor)
 selected=ages[ages.donor_id.isin(chosen)].copy();assert len(selected)==200
 selected.to_csv(O/'G1_prime_selected_200_donors.csv',index=False)
 obs=obs[obs.donor_id.isin(chosen)&obs.broad_type.notna()&obs.is_primary_data&obs.raw_sum.gt(0)].sort_values('soma_joinid')
 obs.merge(selected[['donor_id','age']],on='donor_id',validate='many_to_one').to_parquet(W/'cache/OneK1K_selected_cells.parquet',index=False)
 counts=obs.groupby(['donor_id','broad_type'],observed=True).size().rename('n_cells').reset_index().merge(selected[['donor_id','age']],on='donor_id')
 counts.to_csv(O/'G1_prime_selected_donor_cell_counts.csv',index=False)
 for ct,g in counts.groupby('broad_type'):
  assert g.donor_id.nunique()>=100 and g.age.max()-g.age.min()>=40,(ct,len(g))
 coverage=pd.read_csv(O/'G1_prime_gene_coverage.csv');assert coverage[coverage.dataset.eq('OneK1K')].G1_prime_coverage_pass.all()
 mapping=pd.read_csv(O/'gene_mapping_OneK1K.csv').fillna('');co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');need=set(co.loc[co.coefficient.ne(0)&co.gene.ne('(Intercept)'),'gene'])
 genes=sorted({s for r in mapping[mapping.model_ensembl.isin(need)].itertuples() for s in r.mapped_feature_ids.split('|') if s})
 var=pd.read_csv(O/'OneK1K_analysis_var.csv');v=var[var['_index'].isin(genes)].sort_values('soma_joinid');assert len(v)==len(genes)
 plan={'seed':SEED,'census_release':RELEASE,'dataset_id':'3faad104-2ab8-4434-816d-474d8d2641db','analysis_dataset_version_id':'a3f5651f-cd1a-4d26-8165-74964b79b4f2',
  'n_available_donors':len(ages),'n_selected_donors':len(selected),'age_min':int(selected.age.min()),'age_max':int(selected.age.max()),
  'selection':'40 donors sampled without replacement within each of five observed-age strata; observed age endpoints forced into their original strata before any prediction.',
  'celltype_crosswalk':CROSSWALK,'n_selected_cells':len(obs),'n_requested_unique_features':len(genes),
  'normalization':'log1p(CP10K) using Census raw_sum over the complete released dataset feature universe; source QC retained; ambiguous lineages excluded.',
  'native_small_group_rule':'For 1–15 cells, sample 15 with replacement; for >15 sample without replacement. No outcome-based donor/cell exclusion.',
  'unmeasured_features':'Native zero imputation for prediction; expression summaries remain NA, not biological zero.',
  'pseudocells':100,'pseudocell_size':15,'G1_aggregation':'within-donor median','created_before_prediction_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 (O/'G1_prime_execution_plan.json').write_text(json.dumps(plan,indent=2,ensure_ascii=False));print('PLAN',plan['n_selected_cells'],len(genes),flush=True)
 manifests=[]
 with c.open_soma(census_version=RELEASE) as f:
  for i,lo in enumerate(range(0,len(obs),40000)):
   rr=obs.iloc[lo:lo+40000];path=D/f'OneK1K_{i:03d}.h5ad'
   if not path.exists():
    print('FETCH',i,len(rr),flush=True)
    a=c.get_anndata(f,'Homo sapiens',obs_coords=rr.soma_joinid.to_numpy(dtype=np.int64),var_coords=v.soma_joinid.to_numpy(dtype=np.int64),
     obs_column_names=['soma_joinid','donor_id','cell_type','raw_sum'],var_column_names=['soma_joinid','feature_id','feature_name'])
    assert a.n_obs==len(rr) and a.n_vars==len(v)
    x=a.X;assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
    assert np.all(np.asarray(x.sum(axis=1)).ravel()<=a.obs.raw_sum.to_numpy()+1e-5)
    a.write_h5ad(path.with_suffix('.part.h5ad'),compression='gzip');path.with_suffix('.part.h5ad').rename(path);del a
   check=ad.read_h5ad(path,backed='r');assert check.n_obs==len(rr) and check.n_vars==len(v);check.file.close()
   manifests.append({'path':str(path.relative_to(P)),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'n_cells':len(rr)})
   (O/'OneK1K_count_download_manifest.json').write_text(json.dumps({'plan':plan,'chunks':manifests},indent=2))
   print('DONE',i,path.stat().st_size,flush=True)
 print('ALL_DONE',len(manifests),flush=True)
if __name__=='__main__':main()
