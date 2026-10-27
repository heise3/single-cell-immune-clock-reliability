#!/usr/bin/env python3
"""One local count pass,100 conditional resamples; no new clock fitting."""
from pathlib import Path
from collections import defaultdict
import sys,json,gc,time
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';W=P/'work/round5';A=P/'outputs/round2_A1';R=P/'outputs/round4';sys.path.insert(0,str(A));from run_G1_prime import rng_for,native_indices
from compute_reliability import variance_components,lodo
def noise():
 done=W/'P3_noise_complete.json'
 if done.exists():return pd.read_csv(O/'P3_a_e_conditional_noise.csv')
 man=[]
 for path in [R/'R2_K0_S_evaluation_raw_manifest.json',P/'outputs/round3/R2_E_evaluation_raw_manifest.json']:
  z=json.loads(path.read_text());assert z['complete'];man+=z['chunks']
 models={(ct,v):dict(np.load(R/'models'/f'{ct}_{v}.npz')) for ct in ['B','CD4T','CD8T','MONO','NK'] for v in ['a','e']};score=defaultdict(list);counts=defaultdict(list);libs=defaultdict(list);start=time.monotonic()
 for j,row in enumerate(man):
  a=ad.read_h5ad(P/row['path']);fi={f:i for i,f in enumerate(a.var.feature_id.astype(str))}
  for ct in ['B','CD4T','CD8T','MONO','NK']:
   ix=np.flatnonzero(a.obs.celltype.astype(str).eq(ct));
   if not len(ix):continue
   ma=models[ct,'a'];me=models[ct,'e'];ba=ma['coefficients'];be=me['coefficients'];fa=ma['feature_ids'].astype(str)[ba!=0];fe=me['feature_ids'].astype(str)[be!=0];ba=ba[ba!=0];be=be[be!=0];lib=a.obs.raw_sum.iloc[ix].to_numpy(float);xx=a.X[ix][:,[fi[f] for f in fa]].astype(float).multiply((10000/lib)[:,None]).tocsr();xx.data=np.log1p(xx.data);y=np.asarray(xx@ba).ravel()+float(ma['intercept']);xe=a.X[ix][:,[fi[f] for f in fe]].tocsr();gs=a.obs.gsm.iloc[ix].astype(str).to_numpy()
   for gsm in np.unique(gs):
    ii=np.flatnonzero(gs==gsm);score[gsm,ct].append(y[ii]);counts[gsm,ct].append(xe[ii]);libs[gsm,ct].append(lib[ii])
  print('P3_READ',j+1,len(man),round(time.monotonic()-start,1),flush=True);del a;gc.collect()
 pred=pd.read_csv(R/'F2_a_e_full68_predictions.csv');lk=pred.set_index(['gsm','celltype','variant']).prediction;rows=[];boots=[]
 with threadpool_limits(limits=4):
  for n,(key,ys) in enumerate(sorted(score.items())):
   gsm,ct=key;y=np.concatenate(ys);x=sparse.vstack(counts[key]).tocsr().astype(float);lib=np.concatenate(libs[key]);assert len(y)==x.shape[0]==len(lib);me=models[ct,'e'];b=me['coefficients'];b=b[b!=0];inter=float(me['intercept']);e0=float(np.log2(1+np.asarray(x.sum(axis=0)).ravel()*1e6/lib.sum())@b+inter);assert abs(e0-lk[gsm,ct,'e'])<1e-9,(key,e0,lk[gsm,ct,'e']);rng=rng_for('R5_cell_bootstrap_ae',gsm,ct);multi=rng.multinomial(len(y),np.full(len(y),1/len(y)),size=100);ab=multi@y/len(y);eb=np.log2(1+np.asarray(x.T@multi.T).T*1e6/(multi@lib)[:,None])@b+inter;ix=native_indices(len(y),rng_for('R5_a_algorithm',gsm,ct),rows=10000).reshape(100,100,15);algo=y[ix].mean(axis=(1,2));varpop=float(np.var(y));exact_alg=varpop/1500*((len(y)-15)/(len(y)-1) if len(y)>15 else 1.)
   rows.extend([{'gsm':gsm,'celltype':ct,'model':'a','n_cells':len(y),'conditional_mean':float(y.mean()),'cell_sampling_var':float(ab.var(ddof=1)),'algorithm_var':float(algo.var(ddof=1)),'algorithm_variance_exact_QA':exact_alg,'cell_bootstrap_replicates':100,'algorithm_replicates':100},{'gsm':gsm,'celltype':ct,'model':'e','n_cells':len(y),'conditional_mean':e0,'cell_sampling_var':float(eb.var(ddof=1)),'algorithm_var':0.,'algorithm_variance_exact_QA':0.,'cell_bootstrap_replicates':100,'algorithm_replicates':0}]);boots.extend({'gsm':gsm,'celltype':ct,'replicate':i,'a_cell_conditional':float(aa),'e_cell_conditional':float(ee),'a_fixed_cell_algorithm':float(al)} for i,(aa,ee,al) in enumerate(zip(ab,eb,algo)));del counts[key],libs[key]
   if (n+1)%100==0:print('P3_BOOTSTRAP_GROUP',n+1,len(score),round(time.monotonic()-start,1),flush=True)
 result=pd.DataFrame(rows);assert len(result)==3200;result.to_csv(O/'P3_a_e_conditional_noise.csv',index=False);pd.DataFrame(boots).to_csv(W/'P3_conditional_bootstrap_draws.csv',index=False);done.write_text(json.dumps({'status':'COMPLETE','n_groups':1600,'cells':sum(len(np.concatenate(v)) for v in score.values()),'no_new_training_or_download':True,'native_pseudocell_resampling':'100realizations of100x15; same A1 native_indices implementation','e_pseudocell_variance':'structural0 for deterministic pooled counts, not evidence of no technical noise'},indent=2));return result
def main():
 nv=noise();pred=pd.read_csv(O/'fixed_three_model_predictions.csv');published=pd.read_csv(A/'variance_decomposition.csv');published=published[published.clock.eq('sc-ImmuAging')].copy();published['model']='published';out=published.to_dict('records');pnoise=pd.read_csv(A/'main_predictions.csv');pnoise=pnoise[pnoise.clock.eq('sc-ImmuAging')].rename(columns={'algorithm_var':'algorithm_var'});pnoise['model']='published';noises=pd.concat([nv,pnoise[['gsm','celltype','model','n_cells','conditional_mean','cell_sampling_var','algorithm_var']]],ignore_index=True);noises.to_csv(O/'P3_all_model_noise_inputs.csv',index=False)
 with threadpool_limits(limits=4):
  for (ct,model),g in pred[pred.model.isin(['a','e'])].groupby(['celltype','model']):
   g=g.sort_values(['donor','gsm']).reset_index(drop=True).merge(nv[nv.model.eq(model)],on=['gsm','celltype','model'],validate='one_to_one');acc=lodo(g,g.prediction.to_numpy())
   for metric,y in [('prediction',g.prediction.to_numpy()),('age_adjusted',acc)]:
    z=variance_components(g,y);z.update(celltype=ct,model=model,clock=model,metric=metric);assert z['fit_converged'];out.append(z)
   print('P3_REML',ct,model,flush=True)
 result=pd.DataFrame(out);result.to_csv(O/'P3_three_model_variance_components.csv',index=False);rows=[];components=['donor','visit_nested_plus_unmodelled_residual','batch','cell_sampling','pseudocell_algorithm']
 for (ct,metric),g in result.groupby(['celltype','metric']):
  q=g.set_index('model')
  for comp in components:
   av=q.loc['a',comp+'_variance'];ev=q.loc['e',comp+'_variance'];pv=q.loc['published',comp+'_variance'];rows.append({'celltype':ct,'metric':metric,'component':comp,'published_variance':pv,'a_variance':av,'e_variance':ev,'absolute_reduction_a_minus_e':av-ev,'relative_reduction_a_to_e':1-ev/av if av>0 else np.nan,'published_fraction':q.loc['published',comp+'_fraction'],'a_fraction':q.loc['a',comp+'_fraction'],'e_fraction':q.loc['e',comp+'_fraction']})
 r=pd.DataFrame(rows);r['absolute_reduction_published_minus_e']=r.published_variance-r.e_variance;r['relative_reduction_published_to_e']=1-r.e_variance/r.published_variance;r.to_csv(O/'P3_component_change_plot_data.csv',index=False);dominant=r[r.component.ne('donor')].sort_values('absolute_reduction_a_minus_e',ascending=False).groupby(['celltype','metric'],sort=False).head(1);dominant.to_csv(O/'P3_largest_non_donor_reduction.csv',index=False);(O/'P3_status.json').write_text(json.dumps({'status':'PASS','rows':len(result),'cell_bootstrap':100,'a_algorithm_replicates':100,'e_pseudocell_variance':'structural0','published_components':'exact A1 rows reused, not refitted under a different estimand','same_A1_variance_model':True,'visit_and_unmodelled_technical_residual_separable':False},indent=2));print('P3_COMPLETE',flush=True)
if __name__=='__main__':main()
