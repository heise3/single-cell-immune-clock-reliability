#!/usr/bin/env python3
"""Prepare fixed E/AIDA inputs only; never use their outcomes for model choice."""
from pathlib import Path
import json,time,gc,sys
import numpy as np,pandas as pd,h5py,anndata as ad,cellxgene_census as c
from scipy import sparse
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';W=P/'work/round3';sys.path.insert(0,str(O));from aggregate_training_inputs import Aggregate,sha,CFG
from prepare_training_design import rng_for
CFG={**CFG,'sm.compute_concurrency_level':8,'sm.io_concurrency_level':16}
def prepare():
 ranks=json.loads((O/'R2_training_gene_ranks.json').read_text());union=sorted(set(j for ct in ranks for j in ranks[ct]['full']));ref=pd.read_csv(O/'R2_reference_genes.csv');genes=ref.iloc[union].copy().reset_index().rename(columns={'index':'reference_index'});genes.to_csv(O/'R2_evaluation_gene_universe.csv',index=False)
 split=pd.read_csv(O/'R2_K0_donor_split.csv');eset=set(split.loc[split.split.eq('E'),'donor']);k0=pd.read_csv(P/'outputs/round2_A1/provenance_inputs/K0_primary_with_released_age.csv');k0=k0[k0.donor.isin(eset)];e=pd.read_parquet(P/'work/round2_A1/cache/HIHA_K0_clock_cells.parquet');e=e[e.donor.isin(eset)].rename(columns={'clock_type':'celltype'});e=e.merge(k0[['gsm','sample.subjectAgeAtDraw']].rename(columns={'sample.subjectAgeAtDraw':'age'}),on='gsm',validate='many_to_one');e['split']='E';e['CV_fold']=-1;e['group']=e.gsm+'|'+e.celltype
 source=pd.read_parquet(W/'cache/AIDA_source_obs.parquet');source['age']=pd.to_numeric(source.development_stage.astype(str).str.extract(r'^(\d+)-year-old stage$',expand=False),errors='coerce');don=source.groupby('donor_id').agg(n_age_values=('age','nunique'),age=('age','first'),n_source_cells=('age','size'),development_stage=('development_stage','first'),ethnicity=('self_reported_ethnicity','first')).reset_index();don['commercial_control']=don.donor_id.str.contains('LONZA',case=False,regex=False);don['eligible']=~don.commercial_control&don.age.notna()&don.n_age_values.eq(1);don['reason']=np.select([don.commercial_control,don.age.isna(),don.n_age_values.ne(1)],['COMMERCIAL_COMMON_CONTROL','NO_EXACT_YEAR_AGE','NONUNIQUE_DONOR_AGE'],default='INCLUDE');don.to_csv(O/'R2_AIDA_donor_eligibility.csv',index=False)
 a=pd.read_parquet(W/'cache/AIDA_census_obs.parquet');assert len(a)==len(source) and set(a.donor_id.astype(str))==set(don.donor_id.astype(str));sourcecounts=source.donor_id.astype(str).value_counts();assert (a.donor_id.astype(str).value_counts().sort_index()==sourcecounts.sort_index()).all();cross=json.loads((P/'outputs/round2_A1/G1_prime_execution_plan.json').read_text())['celltype_crosswalk'];cross.update({'CD16-positive, CD56-dim natural killer cell, human':'NK','CD8-positive, alpha-beta cytotoxic T cell':'CD8T','CD8-positive, alpha-beta memory T cell':'CD8T','B cell':'B','mature B cell':'B'});a['celltype']=a.cell_type.astype(str).map(cross);a=a[a.donor_id.isin(don.loc[don.eligible,'donor_id'])&a.celltype.notna()&a.disease.astype(str).eq('normal')&a.raw_sum.gt(0)].copy();a=a.rename(columns={'donor_id':'donor'}).merge(don.loc[don.eligible,['donor_id','age']].rename(columns={'donor_id':'donor'}),on='donor',validate='many_to_one');a['gsm']=a.donor;a['group']=a.donor+'|'+a.celltype;a['split']='AIDA_external';a['CV_fold']=-1
 for name,obs in [('E',e),('AIDA',a)]:obs.sort_values('soma_joinid').to_parquet(W/'cache'/f'{name}_evaluation_cells_pre_QC.parquet',index=False)
 pairs=pd.read_csv(P/'outputs/round2_A1/provenance_inputs/K0_primary_pair_ledger.csv');pairs=pairs[pairs.donor.isin(eset)&pairs.stratum.isin(['75-105d','330-400d'])];pairs.to_csv(O/'R2_E_pair_ledger.csv',index=False);summary={'n_E_donors':e.donor.nunique(),'n_E_visits':e.gsm.nunique(),'E_pairs_by_stratum':pairs.groupby('stratum').size().to_dict(),'AIDA_eligible_noncommercial_donors':int(don.eligible.sum()),'AIDA_commercial_controls_excluded':int(don.commercial_control.sum()),'AIDA_exact_age_min':float(don.loc[don.eligible,'age'].min()),'AIDA_exact_age_max':float(don.loc[don.eligible,'age'].max()),'AIDA_is_primary_data_policy':'Use all cells of fixed v2, including globally nonprimary duplicates of v1; no v1 pooled into analysis. Single v2 observations unique.','AIDA_Census_source_identity':'Exact fixed version; all donor IDs, donor cell counts and total cells match; individual observation_joinid strings not reread from source','E_AIDA_outcomes_used_for_training':False,'n_potential_model_features':len(genes)};(O/'R2_evaluation_input_plan.json').write_text(json.dumps(summary,indent=2));print('EVAL_PLAN',summary,flush=True);return genes,e,a
def run(name,obs,genes,globalvar):
 complete=W/f'{name}_evaluation_input_complete.json'
 if complete.exists():print('REUSE_EVALUATION_INPUT',name,flush=True);return
 obs=obs.sort_values('soma_joinid').reset_index(drop=True);agg=Aggregate(obs,genes.feature_id.to_numpy());vid=genes.soma_joinid.to_numpy(np.int64);mito=globalvar.loc[globalvar.feature_name.str.startswith('MT-',na=False),'soma_joinid'].to_numpy(np.int64);D=W/f'{name}_evaluation_counts';D.mkdir(exist_ok=True);man=[]
 with c.open_soma(census_version='2025-11-08',tiledb_config=CFG) as f:
  X=f['census_data']['homo_sapiens'].ms['RNA'].X['raw']
  for i,lo in enumerate(range(0,len(obs),40000)):
   rr=obs.iloc[lo:lo+40000].copy();file=D/f'{name}_{i:03d}.h5ad';qcfile=file.with_suffix('.QC.json');start=time.monotonic()
   if file.exists() and qcfile.exists():a=ad.read_h5ad(file);x=a.X;o=a.obs.reset_index(drop=True);qc=json.loads(qcfile.read_text())
   else:
    xs=[];os=[];expected=rr.set_index('soma_joinid')
    for mat,(rid,_) in X.read(coords=(rr.soma_joinid.to_numpy(np.int64),slice(None))).blockwise(axis=0,size=5000,reindex_disable_on_axis=1,eager=True).scipy():
     rid=np.asarray(rid,np.int64);o=expected.loc[rid].reset_index();lib=np.asarray(mat.sum(axis=1)).ravel();assert np.array_equal(lib,o.raw_sum.to_numpy());nf=np.diff(mat.tocsr().indptr);mt=np.asarray(mat[:,mito].sum(axis=1)).ravel()/lib;keep=(nf>=200)&(nf<=5000)&(mt<=.25);o=o[keep].copy();o['native_n_features']=nf[keep];o['native_mito_fraction']=mt[keep];os.append(o);xs.append(mat[keep][:,vid].tocsr())
    x=sparse.vstack(xs).tocsr();o=pd.concat(os,ignore_index=True);o.index=o.soma_joinid.astype(str);o.index.name=None;v=genes.copy();v.index=v.feature_id;v.index.name=None;a=ad.AnnData(x,obs=o,var=v);a.write_h5ad(file.with_suffix('.part.h5ad'),compression='lzf');file.with_suffix('.part.h5ad').rename(file);qc={'input_cells':len(rr),'kept_cells':len(o),'excluded_cells':len(rr)-len(o),'full_row_sum_exact':True};qcfile.write_text(json.dumps(qc,indent=2))
   agg.add(x,o);man.append({'path':str(file.relative_to(P)),'bytes':file.stat().st_size,'sha256':sha(file),**qc});(O/f'R2_{name}_evaluation_raw_manifest.json').write_text(json.dumps({'complete':False,'chunks':man},indent=2));print('EVAL_COUNT_CHUNK',name,i,len(o),round(time.monotonic()-start,1),flush=True);del a,x,o;gc.collect()
 agg.save(name);(O/f'R2_{name}_evaluation_raw_manifest.json').write_text(json.dumps({'complete':True,'chunks':man},indent=2));groups=pd.read_csv(W/f'{name}_groups.csv');master=[]
 for row in man:
  a=ad.read_h5ad(P/row['path'],backed='r');master.append(a.obs[['group','soma_joinid','raw_sum']].reset_index(drop=True));a.file.close()
 master=pd.concat(master,ignore_index=True);master.to_parquet(W/'cache'/f'{name}_evaluation_master_obs.parquet',index=False);lookup=master.groupby('group',sort=False).indices;mr=[];mc=[]
 for i,key in enumerate(groups.group):
  cells=np.asarray(lookup[key]);rng=rng_for('R3_evaluation_pseudocell',name,key);draw=np.array([rng.choice(cells,15,replace=len(cells)<=15) for _ in range(100)]);mr.extend(np.repeat(i,1500));mc.extend(draw.ravel())
 G=sparse.csr_matrix((np.full(len(mr),1/1500.),(mr,mc)),shape=(len(groups),len(master)));assert np.allclose(np.asarray(G.sum(axis=1)).ravel(),1);sparse.save_npz(W/f'{name}_native_pseudocell_mean_cell_weights.npz',G);point=np.zeros((len(groups),len(genes)));lo=0
 for i,row in enumerate(man):
  a=ad.read_h5ad(P/row['path']);hi=lo+a.n_obs;x=a.X.astype(np.float64).multiply((10000/a.obs.raw_sum.to_numpy())[:,None]).tocsr();x.data=np.log1p(x.data);z=(G[:,lo:hi]@x).tocoo();point[z.row,z.col]+=z.data;lo=hi;print('EVAL_POINT_CHUNK',name,i,flush=True)
 with h5py.File(W/f'{name}_aggregates.h5','a') as h:
  if 'pseudocell_mean_logCP10K' in h:del h['pseudocell_mean_logCP10K']
  h.create_dataset('pseudocell_mean_logCP10K',data=point,compression='lzf')
 complete.write_text(json.dumps({'status':'COMPLETE','dataset':name,'n_donors':groups.donor.nunique(),'n_visits':groups.gsm.nunique(),'n_cells':len(master),'QC_excluded':len(obs)-len(master),'point_computation':'Exact linear identity: mean of100 predictions on15-cell pseudocells equals prediction of weighted mean expression, using same1500 cell draws; PCA projection also linear','not_used_for_model_selection':True},indent=2));print('EVALUATION_INPUT_COMPLETE',name,len(groups),flush=True)
def main():
 assert (O/'R2_training_gene_ranks.json').exists() and (W/'cache/AIDA_source_obs.parquet').exists();genes,e,a=prepare();v=pd.read_csv(O/'Census_global_var.csv');run('E',e,genes,v);run('AIDA',a,genes,v)
if __name__=='__main__':main()
