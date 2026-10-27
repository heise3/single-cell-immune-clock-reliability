#!/usr/bin/env python3
"""Native linear predictions and fixed-panel RAM, preserving cell-level scores."""
from pathlib import Path
from collections import defaultdict
import json,sys
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
P=Path(__file__).resolve().parents[3];O=P/'outputs/round2_A1/sensitivity87';W=P/'work/round2_A1/sensitivity87';D=W/'main_scores';D.mkdir(exist_ok=True)
sys.path.insert(0,str(O));from run_G1_prime import rng_for,native_indices

def resample_variance(x,rng,B=100):
 return np.array([x[rng.integers(0,len(x),len(x))].mean() for _ in range(B)])
def main():
 population=pd.read_csv(O/'sensitivity87_visit_ledger.csv');allowed_gsms=set(population.gsm)
 manifest=json.loads((O/'main_Census_count_manifest.json').read_text());assert sum(x['n_cells'] for x in manifest['chunks'])==manifest['plan']['n_cells']
 coef=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');mp=pd.read_csv(O/'gene_mapping_HIHA_Census.csv',keep_default_na=False).set_index('model_ensembl');ram=pd.read_csv(O/'RAM_mapped_gene_panels.csv',keep_default_na=False);sc=defaultdict(list);rx=defaultdict(list)
 native=[];checked_ct=set()
 for z in manifest['chunks']:
  a=ad.read_h5ad(P/z['path']);a=a[a.obs.gsm.isin(allowed_gsms)].copy();assert a.n_obs==z['n_cells'];ids=a.var.feature_id.astype(str).tolist();fi={s:i for i,s in enumerate(ids)};lib=a.obs.raw_sum.to_numpy(dtype=float);x=a.X.astype(np.float64).multiply(10000/lib[:,None]).tocsr();x.data=np.log1p(x.data)
  for ct,g in coef[coef.gene.ne('(Intercept)')&coef.coefficient.ne(0)].groupby('celltype'):
   ix=np.flatnonzero(a.obs.clock_type.eq(ct));rr=[];cc=[]
   for j,e in enumerate(g.gene):
    for f in mp.loc[e,'mapped_feature_ids'].split('|'):
     if f:rr.append(fi[f]);cc.append(j)
   A=sparse.csc_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(ids),len(g)));mapped_raw=(a.X[ix]@A).tocsr().astype(np.float64);log=mapped_raw.multiply(10000/lib[ix,None]).tocsr();log.data=np.log1p(log.data);inter=float(coef.loc[coef.celltype.eq(ct)&coef.gene.eq('(Intercept)'),'coefficient'].iloc[0]);y=np.asarray(log@g.coefficient.to_numpy()).ravel()+inter
   for gsm in a.obs.iloc[ix].gsm.unique():
    jj=np.flatnonzero(a.obs.iloc[ix].gsm.to_numpy()==gsm);sc[gsm,ct].append(y[jj])
   if ct not in checked_ct and len(ix):
    for k in range(min(3,len(ix))):native.append({'celltype':ct,'python_prediction':y[k],**dict(zip(g.gene,log[k].toarray().ravel()))})
    checked_ct.add(ct)
  for (ct,direction),panel in ram.groupby(['RAM_celltype','direction'],sort=True):
   ix=np.flatnonzero(a.obs.AIFI_L3.eq(ct));genes=panel.Gene.tolist();rr=[];cc=[]
   for j,s in enumerate(panel.mapped_feature_ids):
    for e in s.split('|'):
     if e:rr.append(fi[e]);cc.append(j)
   A=sparse.csc_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(ids),len(genes)));expr=(x[ix]@A).tocsr()
   # Summing raw counts before log is needed for aliases with >1 feature.
   # RAM mapping has one ID per signature gene in this release; enforce that.
   assert not panel.mapped_feature_ids.str.contains('\\|').any()
   for gsm in a.obs.iloc[ix].gsm.unique():
    jj=np.flatnonzero(a.obs.iloc[ix].gsm.to_numpy()==gsm);rx[gsm,ct,direction].append(expr[jj])
  print('SCORED_CHUNK',z['path'],flush=True)
 # Keep raw cell scores and RAM expression to permit true cell resampling.
 k0=pd.read_csv(O/'sensitivity87_visit_ledger.csv').set_index('gsm');out=[];boots=[]
 for (gsm,ct),arr in sorted(sc.items()):
  y=np.concatenate(arr);r=k0.loc[gsm];np.savez_compressed(D/f'{gsm}_scImmuAging_{ct}.npz',scores=y)
  point=float(y[native_indices(len(y),rng_for('main_point',gsm,ct))].mean());cell=resample_variance(y,rng_for('main_cell_boot',gsm,ct));ix=native_indices(len(y),rng_for('main_algorithm',gsm,ct),rows=10000).reshape(100,100,15);algo=y[ix].mean(axis=(1,2))
  out.append({'clock':'sc-ImmuAging','celltype':ct,'gsm':gsm,'donor':r.donor,'age':r['sample.subjectAgeAtDraw'],'batch':r['batch id'],'n_cells':len(y),'prediction':point,'conditional_mean':y.mean(),'cell_sampling_var':cell.var(ddof=1),'algorithm_var':algo.var(ddof=1)})
  for b,(c,a) in enumerate(zip(cell,algo)):boots.append({'clock':'sc-ImmuAging','celltype':ct,'gsm':gsm,'replicate':b+1,'cell_conditional_mean':c,'fixed_cells_algorithm_prediction':a})
 ramqa=[];cal=[]
 for (ct,direction),panel in ram.groupby(['RAM_celltype','direction'],sort=True):
  keys=sorted([key for key in rx if key[1:]==(ct,direction)]);mats={key[0]:sparse.vstack(rx[key]).tocsr() for key in keys};means=np.vstack([np.asarray(mats[g].mean(axis=0)).ravel() for g in mats]);mu=means.mean(axis=0);sd=means.std(axis=0,ddof=1);valid=(sd>0)&panel.mapped_feature_ids.ne('').to_numpy()
  for j,gene in enumerate(panel.Gene):cal.append({'celltype':ct,'direction':direction,'gene':gene,'reference_mean':mu[j] if panel.mapped_feature_ids.iloc[j] else np.nan,'reference_SD':sd[j] if panel.mapped_feature_ids.iloc[j] else np.nan,'effective_linear_weight_logCP10K':float(1/sd[j]) if valid[j] else 0.0,'used':bool(valid[j]),'measured_in_source':bool(panel.mapped_feature_ids.iloc[j]),'n_visits':len(mats)})
  weights=np.zeros(len(panel));weights[valid]=1/sd[valid];offset=float(mu@weights)
  for gsm,mat in mats.items():
   y=np.asarray(mat@weights).ravel()-offset;r=k0.loc[gsm];cell=resample_variance(y,rng_for('RAM_cell_boot',gsm,ct,direction));clock='RAM_'+direction;pred=float(y.mean());expected=float(((means[list(mats).index(gsm),valid]-mu[valid])/sd[valid]).sum());assert abs(pred-expected)<1e-7
   sparse.save_npz(D/f'{gsm}_{clock}_{ct.replace("/","_")}.npz',mat);out.append({'clock':clock,'celltype':ct,'gsm':gsm,'donor':r.donor,'age':r['sample.subjectAgeAtDraw'],'batch':r['batch id'],'n_cells':len(y),'prediction':pred,'conditional_mean':pred,'cell_sampling_var':cell.var(ddof=1),'algorithm_var':0.0})
   for b,c in enumerate(cell):boots.append({'clock':clock,'celltype':ct,'gsm':gsm,'replicate':b+1,'cell_conditional_mean':c,'fixed_cells_algorithm_prediction':pred})
  ramqa.append({'celltype':ct,'direction':direction,'n_visits':len(mats),'n_selected_genes':len(panel),'n_used_genes':int(valid.sum()),'reference_score_intercept':-offset,'undefined_or_absent_genes':'|'.join(panel.loc[~valid,'Gene']),'pseudocell_variance_status':'NOT_APPLICABLE_NATIVE_DETERMINISTIC'})
 result=pd.DataFrame(out);assert set(result.gsm).issubset(allowed_gsms);assert set(result.loc[result.clock.eq('sc-ImmuAging'),'gsm'])==allowed_gsms;assert result[result.clock.eq('sc-ImmuAging')].donor.nunique()==population.donor.nunique();result.groupby(['clock','celltype']).agg(n_visits=('gsm','nunique'),n_donors=('donor','nunique'),cells_min=('n_cells','min'),cells_median=('n_cells','median')).to_csv(O/'clock_sample_coverage.csv');result.to_csv(O/'main_predictions.csv',index=False);pd.DataFrame(boots).to_csv(O/'cell_and_algorithm_bootstrap.csv',index=False);pd.DataFrame(cal).to_csv(O/'RAM_reference_calibration.csv',index=False);pd.DataFrame(ramqa).to_csv(O/'RAM_runtime_gene_QA.csv',index=False);pd.DataFrame(native).to_csv(O/'main_native_prediction_input.csv',index=False)
 print('ALL_CLOCKS_SCORED',len(out),flush=True)
if __name__=='__main__':main()
