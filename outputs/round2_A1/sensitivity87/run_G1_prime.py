#!/usr/bin/env python3
"""A1 donor-level matched pseudocell masking; no age-range tuning after prediction."""
from pathlib import Path
from collections import defaultdict
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,stats
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';SEED=20261025
def rng_for(*args):
 h=int.from_bytes(hashlib.sha256('|'.join(map(str,args)).encode()).digest()[:8],'little')
 return np.random.default_rng(np.random.SeedSequence([SEED,h]))
def native_indices(n,rng,rows=100,size=15):
 # Sequential rejection gives independent uniform SRSWOR rows, without the
 # incorrect global-replace=False behavior of choice(size=(rows,size)).
 z=rng.integers(0,n,size=(rows,size))
 if n<=size:return z
 for j in range(1,size):
  bad=np.any(z[:,:j]==z[:,j,None],axis=1)
  while bad.any():
   z[bad,j]=rng.integers(0,n,size=int(bad.sum()));bad=np.any(z[:,:j]==z[:,j,None],axis=1)
 return z
def main():
 plan=json.loads((O/'G1_prime_execution_plan.json').read_text());manifest=json.loads((O/'OneK1K_count_download_manifest.json').read_text());paths=[P/r['path'] for r in manifest['chunks']]
 assert sum(r['n_cells'] for r in manifest['chunks'])==plan['n_selected_cells'],'Incomplete download'
 co=pd.read_csv(P/'outputs/scImmuAging_coefficients.csv');mapping=pd.read_csv(O/'gene_mapping_OneK1K.csv').fillna('');gm=mapping.set_index('model_ensembl')
 residual=pd.read_csv(O/'G1_prime_residual_genes_pending_expression.csv').fillna('');ages=pd.read_csv(O/'G1_prime_selected_200_donors.csv').set_index('donor_id').age.to_dict()
 cross=plan['celltype_crosswalk'];ct_models={};scores=defaultdict(list);expr_sum=defaultdict(lambda:None);cp_sum=defaultdict(lambda:None);ncell=defaultdict(int);native_rows=[]
 for ct,g in co[co.gene.ne('(Intercept)')&co.coefficient.ne(0)].groupby('celltype'):
  genes=g.gene.tolist();ct_models[ct]={'genes':genes,'beta':g.coefficient.to_numpy(),'intercept':float(co.loc[co.celltype.eq(ct)&co.gene.eq('(Intercept)'),'coefficient'].iloc[0]),'missing':set(residual.loc[residual.celltype.eq(ct),'gene'])}
 for ip,path in enumerate(paths):
  a=ad.read_h5ad(path);feature_ids=a.var.feature_id.astype(str).to_list();fi={s:i for i,s in enumerate(feature_ids)};types=a.obs.cell_type.astype(str).map(cross);lib=a.obs.raw_sum.to_numpy(dtype=float)
  for ct,model in ct_models.items():
   ii=np.flatnonzero(types.eq(ct));
   if not len(ii):continue
   genes=model['genes'];rr=[];cc=[]
   for j,e in enumerate(genes):
    for f in gm.loc[e,'mapped_feature_ids'].split('|'):
     if f:rr.append(fi[f]);cc.append(j)
   A=sparse.csc_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(fi),len(genes)))
   counts=(a.X[ii]@A).tocsr().astype(np.float64);cp=counts.multiply(10000/lib[ii,None]).tocsr();log=cp.copy();log.data=np.log1p(log.data)
   miss=np.array([e in model['missing'] for e in genes]);full=np.asarray(log@model['beta']).ravel()+model['intercept'];delta=np.asarray(log[:,miss]@model['beta'][miss]).ravel()
   donors=a.obs.iloc[ii].donor_id.astype(str).to_numpy();mg=[j for j,e in enumerate(genes) if miss[j]]
   for donor in np.unique(donors):
    jj=np.flatnonzero(donors==donor);key=(donor,ct);scores[key].append(np.c_[full[jj],delta[jj]])
    sums=np.asarray(log[jj][:,mg].sum(axis=0)).ravel();cps=np.asarray(cp[jj][:,mg].sum(axis=0)).ravel()
    expr_sum[key]=sums if expr_sum[key] is None else expr_sum[key]+sums;cp_sum[key]=cps if cp_sum[key] is None else cp_sum[key]+cps;ncell[key]+=len(jj)
   # Export real normalized feature vectors to check native glmnet prediction.
   if ip==0:
    z=log[:3].toarray()
    for k in range(len(z)):
     native_rows.append({'celltype':ct,'python_prediction':float(full[k]),**dict(zip(genes,z[k]))})
  print('SCORED',ip,a.n_obs,flush=True);del a
 pd.DataFrame(native_rows).to_csv(O/'G1_prime_native_glmnet_validation_input.csv',index=False)
 donor_rows=[];gene_rows=[]
 for key,arrs in sorted(scores.items()):
  donor,ct=key;x=np.concatenate(arrs);ix=native_indices(len(x),rng_for('G1_prime',donor,ct));ps=x[ix].mean(axis=1)
  f=float(np.median(ps[:,0]));m=float(np.median(ps[:,0]-ps[:,1]));donor_rows.append({'donor':donor,'celltype':ct,'age':ages[donor],'n_cells':len(x),'full_median':f,'masked_median':m,'full_conditional_mean':float(x[:,0].mean()),'missing_conditional_mean_component':float(x[:,1].mean()),'native_small_group_replacement':len(x)<=15})
  genes=[e for e in ct_models[ct]['genes'] if e in ct_models[ct]['missing']]
  for j,e in enumerate(genes):
   measured=bool(gm.loc[e,'mapped_feature_ids']);gene_rows.append({'donor':donor,'celltype':ct,'gene':e,'n_cells':len(x),'mean_log1p_CP10K':expr_sum[key][j]/ncell[key] if measured else np.nan,'mean_CP10K':cp_sum[key][j]/ncell[key] if measured else np.nan,'measured':measured})
 donors=pd.DataFrame(donor_rows);donors.to_csv(O/'G1_prime_donor_predictions.csv',index=False);gx=pd.DataFrame(gene_rows);gx.to_csv(O/'G1_prime_missing_gene_donor_expression.csv',index=False)
 summary=gx.groupby(['celltype','gene'],observed=True).agg(expression_mean_log1p_CP10K=('mean_log1p_CP10K','mean'),between_donor_SD_log1p_CP10K=('mean_log1p_CP10K','std'),expression_mean_CP10K=('mean_CP10K','mean'),between_donor_SD_CP10K=('mean_CP10K','std'),n_measured_donors=('mean_log1p_CP10K','count')).reset_index()
 metadata=residual.drop(columns=['expression_mean','expression_between_donor_SD','expression_status'],errors='ignore');metadata.merge(summary,on=['celltype','gene'],how='left').to_csv(O/'G1_prime_missing_gene_list.csv',index=False)
 out=[];pre=[]
 for ct,g in donors.groupby('celltype'):
  y=g.full_median.to_numpy();z=g.masked_median.to_numpy();d=y-z;u=g.full_conditional_mean.to_numpy();w=g.missing_conditional_mean_component.to_numpy();vu=np.var(u,ddof=1);vy=np.var(y,ddof=1);vz=np.var(z,ddof=1)
  pre.append({'celltype':ct,'n_donors':len(g),'conditional_mean_full_variance':vu,'missing_component_variance':np.var(w,ddof=1),'missing_component_variance_fraction':np.var(w,ddof=1)/vu,'missing_covariance_allocated_fraction':np.cov(u,w,ddof=1)[0,1]/vu,'variance_removed_fraction':(vu-np.var(u-w,ddof=1))/vu,'note':'Exact linear donor-mean covariance decomposition; may be negative or exceed one. Unmeasured genes imputed zero by native prediction, expression not asserted zero.'})
  r=stats.pearsonr(y,z).statistic;rho=stats.spearmanr(y,z).statistic;ratio=vz/vy
  eligible=len(g)>=100 and g.age.max()-g.age.min()>=40
  passed=eligible and r>=.9 and .8<=ratio<=1.25
  out.append({'celltype':ct,'n_donors':len(g),'age_min':g.age.min(),'age_max':g.age.max(),'Pearson_r':r,'Spearman_rho':rho,'variance_ratio_masked_full':ratio,'bias_masked_minus_full_years':float(np.mean(z-y)),'full_between_donor_variance':vy,'masked_between_donor_variance':vz,'median_effective_missing_component_var_fraction':np.var(d,ddof=1)/vy,'median_effective_missing_covariance_fraction':np.cov(y,d,ddof=1)[0,1]/vy,'G1_prime_status':'PASS' if passed else ('EXCLUDE_FROM_MAIN' if eligible else 'NOT_EVALUABLE')})
 pd.DataFrame(pre).to_csv(O/'G1_prime_variance_precheck.csv',index=False);pd.DataFrame(out).to_csv(O/'G1_prime_results.csv',index=False)
 print(pd.DataFrame(out).to_string(index=False),flush=True)
if __name__=='__main__':main()
