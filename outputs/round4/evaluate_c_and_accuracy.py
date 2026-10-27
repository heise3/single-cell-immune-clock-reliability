#!/usr/bin/env python3
"""Secondary repeated OOF predictions and a/c/e accuracy/bias decomposition."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd,h5py
from scipy import stats
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';W=P/'work/round4';R=P/'outputs/round3';V=P/'work/round3';A=P/'outputs/round2_A1';sys.path.insert(0,str(A));from compute_reliability import lodo,pair_array,icc_sem
def accuracy(pred):
 rows=[]
 for (ds,ct,var),g in pred.groupby(['dataset','celltype','variant']):
  g=g.copy();g['error']=g.prediction-g.age;dd=g.groupby('donor').agg(pred=('prediction','mean'),age=('age','mean'),bias=('error','mean'));bias=float(dd.bias.mean());ns=g.groupby('donor').size();w=1/len(dd)/g.donor.map(ns).to_numpy();error=g.error.to_numpy();mse=float(w@(error*error));centered=float(w@((error-bias)**2));mae=float(w@abs(error));frac=bias*bias/mse if mse else np.nan;r=float(stats.pearsonr(dd.pred,dd.age).statistic) if dd.pred.std()>1e-6 else np.nan;rows.append({'dataset':ds,'celltype':ct,'variant':var,'n_donors':len(dd),'n_visits':len(g),'MAE':mae,'Pearson_r':r,'bias_mean_error':bias,'RMSE':np.sqrt(mse),'centered_RMSE_scatter':np.sqrt(centered),'MSE':mse,'bias_squared':bias*bias,'centered_MSE':centered,'bias_squared_fraction_MSE':frac,'diagnostic_MAE_after_removing_global_mean_error':float(w@abs(error-bias)),'dominant_AIDA_error_component':'fixed_offset' if frac>.5 else 'scatter','weighting':'equal donor weight, within-donor visit weight1/n','offset_removal_scope':'diagnostic using external outcomes; not deployed recalibration'});assert np.isclose(mse,bias*bias+centered,atol=1e-9)
 return pd.DataFrame(rows)
def main():
 assert (W/'c_models/all_complete.json').exists() and (O/'F2_decision.json').exists();meta=pd.read_csv(W/'K0_full_groups.csv');am=pd.read_csv(V/'AIDA_groups.csv');genes=pd.read_csv(O/'K0_prediction_gene_universe.csv');lk={f:i for i,f in enumerate(genes.feature_id)};splits=pd.read_csv(O/'c_50x2_donor_fold_ledger.csv');pairs=pd.read_csv(A/'provenance_inputs/K0_primary_pair_ledger.csv');pairs=pairs[pairs.stratum.isin(['75-105d','330-400d'])];repeats=[];ensembles=[];crows=[];oofrows=[]
 with h5py.File(W/'K0_full_inputs.h5') as h,h5py.File(V/'AIDA_aggregates.h5') as ha,threadpool_limits(limits=4):
  full=h['pseudocell_mean_logCP10K'][:];ax=ha['pseudocell_mean_logCP10K'][:]
  for ct in ['B','CD4T','CD8T','MONO','NK']:
   take=meta.celltype.eq(ct).to_numpy();g=meta[take].reset_index(drop=True);x=full[take];at=am.celltype.eq(ct).to_numpy();agm=am[at].reset_index(drop=True);xa=ax[at];one=dict(np.load(W/'OneK_c_stats'/f'{ct}_full.npz'));ol={int(k):i for i,k in enumerate(one['reference_index'])};pk=np.zeros((50,len(g)));pa=np.zeros((50,len(agm)));po=np.zeros((50,len(one['test_donor'])))
   for rep in range(50):
    fold=splits[splits.repeat.eq(rep)].set_index('donor').heldout_fold;ff=g.donor.map(fold).to_numpy();assert np.isfinite(ff).all()
    for held in [0,1]:
     m=np.load(W/'c_models'/f'{ct}_r{rep:02d}_h{held}.npz');b=m['coefficients'];inter=float(m['intercept']);ix=[lk[f] for f in m['feature_ids'].astype(str)];oi=[ol[int(j)] for j in m['reference_index']];test=ff==held;pk[rep,test]=x[test][:,ix]@b+inter;pa[rep]+=(xa[:,ix]@b+inter)/2;po[rep]+=(one['test_mean'][:,oi]@b+inter)/2
    acc=lodo(g,pk[rep]);deg=np.std(acc,ddof=1)<1e-6
    for st,pp in pairs.groupby('stratum'):
     q,a=pair_array(g,acc,pp);icc=icc_sem(a)[0] if not deg else np.nan;repeats.append({'celltype':ct,'repeat':rep,'stratum':st,'n_full_donors':68,'n_pairs':len(q),'ICC2_1':icc,'status':'ESTIMABLE' if np.isfinite(icc) else 'NOT_ESTIMABLE','metric':'LODO_age_acceleration','uncertainty_scope':'one repeated2fold OOF prediction; no donor bootstrap for this repeat'})
    oofrows.extend({'celltype':ct,'repeat':rep,'donor':d,'gsm':s,'age':float(age),'prediction':float(p),'age_acceleration':float(a),'heldout_fold':int(f)} for d,s,age,p,a,f in zip(g.donor,g.gsm,g.age,pk[rep],acc,ff))
   avg=pk.mean(axis=0);res=lodo(g,avg)
   for st,pp in pairs.groupby('stratum'):
    q,a=pair_array(g,res,pp);ensembles.append({'celltype':ct,'stratum':st,'n_donors':68,'n_pairs':len(q),'OOF_ensemble_ICC2_1':icc_sem(a)[0],'scope':'Average50 OOF predictions per donor/visit, then refitfull68LODO; distinct from median50 repeat ICCs'})
   for ds,gm,pred in [('Allen_K0_full',g,avg),('AIDA_external',agm,pa.mean(axis=0))]:crows.extend({'dataset':ds,'celltype':ct,'variant':'c','donor':d,'gsm':s,'age':float(a),'prediction':float(p)} for d,s,a,p in zip(gm.donor,gm.gsm,gm.age,pred))
   crows.extend({'dataset':'OneK1K_test20','celltype':ct,'variant':'c','donor':str(d),'gsm':str(d),'age':float(a),'prediction':float(p)} for d,a,p in zip(one['test_donor'],one['test_age'],po.mean(axis=0)));np.savez_compressed(W/f'c_{ct}_repeat_predictions.npz',Allen_OOF=pk,AIDA_two_fold_ensemble=pa,OneK_two_fold_ensemble=po);print('C_OOF_EVALUATED',ct,flush=True)
 rep=pd.DataFrame(repeats);rep.to_csv(O/'Figure4_c_50_repeat_ICC.csv',index=False);summary=[];base=pd.read_csv(O/'Figure4_F2_full68_a_e_ICC.csv')
 for (ct,st),q in rep.groupby(['celltype','stratum']):
  z=q.ICC2_1.to_numpy();valid=np.isfinite(z);lo,median,hi=np.nanpercentile(z,[2.5,50,97.5]);b=base[base.celltype.eq(ct)&base.variant.eq('a')&base.stratum.eq(st)].ICC2_1.iloc[0];summary.append({'celltype':ct,'stratum':st,'n_repeats':50,'n_estimable_repeats':int(valid.sum()),'ICC_median':median,'ICC_repeat_p025':lo,'ICC_repeat_p975':hi,'baseline_a_full68_ICC':b,'delta_ICC_median':median-b,'delta_repeat_p025':lo-b,'delta_repeat_p975':hi-b,'interval_interpretation':'Empirical50 repeated-split percentiles, NOT a donor-bootstrap95%CI'})
 pd.DataFrame(summary).to_csv(O/'Figure4_c_repeated_crossfit_summary.csv',index=False);pd.DataFrame(ensembles).to_csv(O/'c_OOF_ensemble_ICC.csv',index=False);pd.DataFrame(oofrows).to_csv(O/'c_Allen_50_repeat_OOF_predictions.csv',index=False);c=pd.DataFrame(crows);c.to_csv(O/'c_ensemble_predictions.csv',index=False);old=pd.read_csv(R/'R2_all_predictions.csv');old=old[old.variant.isin(['a','e'])&~old.dataset.eq('Allen_E')];primary=pd.read_csv(O/'F2_a_e_full68_predictions.csv');allpred=pd.concat([old,primary,c],ignore_index=True);allpred.to_csv(O/'a_c_e_all_accuracy_predictions.csv',index=False);a=accuracy(allpred);a.to_csv(O/'Figure5_accuracy_bias_decomposition.csv',index=False);aa=a[a.dataset.eq('AIDA_external')];pool=aa.groupby('variant')[['MSE','bias_squared','centered_MSE']].mean().reset_index();pool['celltype_specific_offset_fraction_MSE']=pool.bias_squared/pool.MSE;pool['interpretation']='Equal-weight mean acrossfive cell types; keeps type-specific offsets rather than cancelling opposite-signed bias';pool.to_csv(O/'AIDA_equal_type_error_summary.csv',index=False);(O/'c_and_accuracy_status.json').write_text(json.dumps({'status':'PASS','c_fitted_models':500,'repeats':50,'folds':2,'c_secondary_only':True,'accuracy_rows':len(a),'c_accuracy_definition':'Allen average50 trueOOF predictions; external average100 opposite-fold-trained models','F2_not_changed_by_c_or_AIDA':True,'percentile_scope':'Repeated-split distribution, not a confidence interval','AIDA_offset_MSE_identity_verified':True},indent=2));print('C_ACCURACY_COMPLETE',flush=True)
if __name__=='__main__':main()
