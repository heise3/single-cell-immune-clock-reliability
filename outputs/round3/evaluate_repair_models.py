#!/usr/bin/env python3
"""Frozen-model E reliability, independent accuracy, and unchanged F decision."""
from pathlib import Path
import sys,json,datetime
import numpy as np,pandas as pd,h5py
from scipy import stats
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';W=P/'work/round3';M=W/'models';SEED=20261025;TYPES=['B','CD4T','CD8T','MONO','NK'];VARS=['a','b','c','d','e'];sys.path.insert(0,str(P/'outputs/round2_A1'));from compute_reliability import lodo,batched_lodo,icc_sem,pair_array
def accuracy(pred):
 out=[]
 for (dataset,ct,var),g in pred.groupby(['dataset','celltype','variant']):
  g=g.copy();g['abs_error']=abs(g.prediction-g.age);dd=g.groupby('donor').agg(prediction=('prediction','mean'),age=('age','mean'),MAE=('abs_error','mean'));r=float(stats.pearsonr(dd.prediction,dd.age).statistic) if dd.prediction.std()>1e-6 and dd.age.std()>0 else np.nan;out.append({'dataset':dataset,'celltype':ct,'variant':var,'n_donors':len(dd),'n_observations':len(g),'MAE':dd.MAE.mean(),'Pearson_r':r,'median_absolute_error':g.groupby('donor').abs_error.median().median(),'age_min':g.age.min(),'age_max':g.age.max(),'MAE_definition':'equal donor weight; mean within-donor visit absolute error','Pearson_definition':'one mean predicted/chronological age per donor','prediction_SD_between_donor_means':dd.prediction.std()})
 return pd.DataFrame(out)
def predict_inputs():
 rows=[];coverage=[];genes=pd.read_csv(O/'R2_evaluation_gene_universe.csv');lookup={f:i for i,f in enumerate(genes.feature_id)};a_var=pd.read_csv(O/'AIDA_source_var.csv');acol='_index' if '_index' in a_var else 'feature_id';apresent=set(a_var[acol].astype(str));one=pd.read_csv(O/'R2_OneK1K_holdout_predictions.csv');one['gsm']=one.donor
 for name in ['E','AIDA']:
  assert (W/f'{name}_evaluation_input_complete.json').exists();meta=pd.read_csv(W/f'{name}_groups.csv')
  with h5py.File(W/f'{name}_aggregates.h5') as h:
   ps=h['pseudocell_mean_logCP10K'][:];pb=h['pseudobulk_logCPM'][:]
   for ct in TYPES:
    mask=meta.celltype.eq(ct).to_numpy();g=meta[mask]
    for var in VARS:
     model=np.load(M/f'{ct}_{var}.npz');fids=model['feature_ids'].astype(str);beta=model['coefficients'];ix=np.array([lookup[f] for f in fids],int);x=(pb if var=='e' else ps)[mask][:,ix];pr=x@beta+float(model['intercept']);rows.extend({'dataset':'Allen_E' if name=='E' else 'AIDA_external','celltype':ct,'variant':var,'donor':d,'gsm':gsm,'age':float(age),'prediction':float(p)} for d,gsm,age,p in zip(g.donor,g.gsm,g.age,pr));nonzero=fids[beta!=0];missing=[f for f in nonzero if f not in apresent] if name=='AIDA' else [];coverage.append({'dataset':name,'celltype':ct,'variant':var,'n_nonzero_features':len(nonzero),'n_missing_nonzero_features':len(missing),'nonzero_coverage':1-len(missing)/len(nonzero) if len(nonzero) else np.nan,'missing_feature_ids':'|'.join(missing),'imputation':'native zero input for annotation-absent feature; not a claim expression is biologically zero'})
 pred=pd.concat([one,pd.DataFrame(rows)],ignore_index=True);pred.to_csv(O/'R2_all_predictions.csv',index=False);pd.DataFrame(coverage).to_csv(O/'R2_evaluation_feature_coverage.csv',index=False);return pred
def reliability(pred):
 pairs=pd.read_csv(O/'R2_E_pair_ledger.csv');e=pred[pred.dataset.eq('Allen_E')];ids=sorted(e.donor.unique());boot=np.random.default_rng(np.random.SeedSequence([SEED,3202])).multinomial(len(ids),np.full(len(ids),1/len(ids)),size=2000);out=[];bs={};point_res=[]
 for (ct,var),g in e.groupby(['celltype','variant'],sort=True):
  g=g.sort_values(['donor','gsm']).reset_index(drop=True);assert sorted(g.donor.unique())==ids;raw=g.prediction.to_numpy();acc=lodo(g,raw);samples=batched_lodo(g,raw,boot);g['age_acceleration']=acc;point_res.append(g);degenerate=np.std(acc,ddof=1)<1e-6
  for st,pp in pairs.groupby('stratum',sort=True):
   q,a=pair_array(g,acc,pp);pt=icc_sem(a)[0] if not degenerate and np.std(a)>1e-6 else np.nan;vals=[]
   for k,m in enumerate(boot):
    qq,z=pair_array(g,samples[k],pp);repeat=np.array([m[ids.index(d)] for d in qq.donor]);z=np.repeat(z,repeat,axis=0);vals.append(icc_sem(z)[0] if np.std(z)>1e-6 else np.nan)
   vals=np.array(vals);finite=np.isfinite(vals);lo,hi=np.nanpercentile(vals,[2.5,97.5]) if finite.any() else [np.nan,np.nan];out.append({'celltype':ct,'variant':var,'stratum':st,'n_E_donors':len(ids),'n_paired_donors':len(q),'metric':'LODO_age_acceleration','ICC2_1':pt,'ICC_CI_low':lo,'ICC_CI_high':hi,'bootstrap_replicates':2000,'bootstrap_valid':int(finite.sum()),'status':'ESTIMABLE' if np.isfinite(pt) else 'NOT_ESTIMABLE_CONSTANT_OR_DEGENERATE_SCORE','CI_scope':'conditional on fixed trained model; refit LODO and resample E donor clusters, same draws all variants'});bs[ct,var,st]=vals
  print('E_RELIABILITY',ct,var,flush=True)
 result=pd.DataFrame(out);deltas=[]
 for r in result.itertuples():
  b=result[(result.celltype==r.celltype)&(result.variant=='a')&(result.stratum==r.stratum)].iloc[0];bv=bs[r.celltype,r.variant,r.stratum]-bs[r.celltype,'a',r.stratum];finite=np.isfinite(bv);lo,hi=np.nanpercentile(bv,[2.5,97.5]) if finite.any() else [np.nan,np.nan];deltas.append({'celltype':r.celltype,'variant':r.variant,'stratum':r.stratum,'delta_ICC_vs_replicated_baseline':r.ICC2_1-b.ICC2_1,'delta_ICC_CI_low':lo,'delta_ICC_CI_high':hi,'paired_bootstrap_valid':int(finite.sum())})
 result=result.merge(pd.DataFrame(deltas),on=['celltype','variant','stratum'],validate='one_to_one');result.to_csv(O/'R2_E_reliability_results.csv',index=False);pd.concat(point_res,ignore_index=True).to_csv(O/'R2_E_predictions_with_LODO.csv',index=False);np.savez_compressed(W/'R2_ICC_bootstrap_distributions.npz',**{'__'.join(k):v for k,v in bs.items()});return result
def main():
 assert (M/'all_models_complete.json').exists();pred=predict_inputs();acc=accuracy(pred);acc.to_csv(O/'R2_accuracy_results.csv',index=False);rel=reliability(pred);one=acc[acc.dataset.eq('OneK1K_test20')][['celltype','variant','MAE']];base=one[one.variant.eq('a')].set_index('celltype').MAE;gate=rel[rel.stratum.eq('75-105d')&rel.celltype.isin(['B','MONO','NK'])].merge(one,on=['celltype','variant'],validate='one_to_one');gate['baseline_OneK1K_MAE']=gate.celltype.map(base);gate['MAE_increase']=gate.MAE-gate.baseline_OneK1K_MAE;gate['ICC_gain_pass']=gate.delta_ICC_vs_replicated_baseline.ge(.15);gate['accuracy_guardrail_pass']=gate.MAE_increase.le(1.5);gate['type_meets_F']=gate.ICC_gain_pass&gate.accuracy_guardrail_pass&gate.variant.ne('a');gate.to_csv(O/'R2_F_gate_by_type.csv',index=False);counts=gate[gate.variant.ne('a')].groupby('variant').type_meets_F.sum();success=counts[counts.ge(2)].index.tolist();verdict='PASS' if success else 'FAIL';summary={'seed':SEED,'repair_verdict':verdict,'successful_variants':success,'n_target_types_meeting_both_rules_by_variant':counts.to_dict(),'rule':'same variant,at least2 of B/MONO/NK90d residual ICC gain>=0.15 AND OneK held-out20% MAE increase<=1.5 years in those types','manuscript_route':'audit plus prespecified repair evidence' if success else 'pure audit; no repair efficacy claim','F_unchanged':True,'E_donors':pred[pred.dataset.eq('Allen_E')].donor.nunique(),'prior_A1_E_audited':'E was previously in A1 descriptive audit; not a pristine new confirmatory cohort. S/E separation prevents this-round feature/fitting leakage.','constant_score_rule':'Residual SD<1e-6 years or degenerate paired variance yields undefined ICC; cannot count as improvement. Not an access failure.','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(O/'R2_F_decision.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False));print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
