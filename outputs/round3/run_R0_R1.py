#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,datetime,warnings
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';A=P/'outputs/round2_A1';SEED=20261025
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def ep(obs,a):
 a=np.asarray(a);a=a[np.isfinite(a)];return min(1.,2*min((np.sum(a<=obs)+1)/(len(a)+1),(np.sum(a>=obs)+1)/(len(a)+1)))
def holm(p):
 p=np.asarray(p);ix=np.argsort(p);q=np.maximum.accumulate(p[ix]*np.arange(len(p),0,-1));out=np.empty(len(p));out[ix]=np.minimum(q,1);return out
def main():
 ledger=json.loads((P/'outputs/gate_ledger.json').read_text());assert ledger['round3_F']['written_before_round3_calculation'];cp=P/'outputs/round2_A1_checkpoint.zip';expected=json.loads((P/'outputs/round2_A1_checkpoint_receipt.json').read_text())['sha256'];assert sha(cp)==expected
 files=pd.read_csv(A/'SHA256SUMS.csv');assert all(sha(A/r.file)==r.sha256 for r in files.itertuples());(O/'checkpoint_recovery_QA.json').write_text(json.dumps({'archive_sha256':expected,'verified_files':len(files),'status':'PASS','preserved_original_files':True},indent=2))
 r=pd.read_csv(A/'reliability_results.csv');x=r[r.clock.eq('sc-ImmuAging')];cols=['ICC2_1','ICC_CI_low','ICC_CI_high','SEM','SEM_CI_low','SEM_CI_high','MDC95','MDC95_CI_low','MDC95_CI_high'];wide=x.pivot(index=['celltype','stratum','n_donors'],columns='metric',values=cols);wide.columns=[f'{metric}__{stat}' for stat,metric in wide.columns];wide=wide.reset_index();wide.to_csv(O/'Table1_scImmuAging_reliability_wide.csv',index=False)
 v=pd.read_csv(A/'variance_decomposition.csv');comp=['donor','visit_nested_plus_unmodelled_residual','batch','cell_sampling','pseudocell_algorithm'];ff=[k+'_fraction' for k in comp];v['largest_component']=v[ff].idxmax(axis=1).str.replace('_fraction','',regex=False);v['largest_component_fraction']=v[ff].max(axis=1);v[['clock','celltype','metric','n_donors','n_visits',*ff,'largest_component','largest_component_fraction']].to_csv(O/'Table2_variance_components_42.csv',index=False)
 counts=v.groupby(['clock','metric','largest_component']).size().rename('n_rows').reset_index();counts.to_csv(O/'R0_dominant_component_counts.csv',index=False);summary={'status':'PASS','Table1_rows':len(wide),'Table2_rows':len(v),'dominance_counts_all42':v.largest_component.value_counts().to_dict(),'scImmuAging_age_adjusted_dominance':v[v.clock.eq('sc-ImmuAging')&v.metric.eq('age_adjusted')].largest_component.value_counts().to_dict(),'interpretation':'Largest individual component is not necessarily more than half of total variance; visit component includes unmodelled residual.'};(O/'R0_summary.json').write_text(json.dumps(summary,indent=2));print('R0',summary,flush=True)
 paired=pd.read_csv(A/'BCG_38_paired_CD8.csv');delta=paired.change.to_numpy();mu=float(delta.mean());nobs=int((delta<0).sum());mobs=float(np.median(delta[delta<0]));B=100000;rng=np.random.default_rng(np.random.SeedSequence([SEED,3101]));z=rng.standard_normal((B,38));out=[]
 for sem in np.round(np.arange(.5,6.00001,.1),1):
  a=mu+np.sqrt(2)*sem*z;n=(a<0).sum(axis=1)
  with warnings.catch_warnings():warnings.filterwarnings('ignore',message='All-NaN slice encountered');med=np.nanmedian(np.where(a<0,a,np.nan),axis=1)
  p1=ep(nobs,n);p2=ep(mobs,med);q=holm([p1,p2]);out.append({'SEM':sem,'null_mean':mu,'null_SD':np.sqrt(2)*sem,'n_simulations':B,'n_nonempty_negative_subgroups':int(np.isfinite(med).sum()),'observed_n_negative':nobs,'observed_negative_subgroup_median':mobs,'negative_count_empirical_p':p1,'negative_count_Holm_p':q[0],'negative_median_empirical_p':p2,'negative_median_Holm_p':q[1],'any_Holm_p_lt_0_05':bool((q<.05).any()),'any_unadjusted_p_lt_0_05':bool(min(p1,p2)<.05),'null_count_p025':np.quantile(n,.025),'null_count_p975':np.quantile(n,.975),'null_negative_median_p025':np.nanquantile(med,.025),'null_negative_median_p975':np.nanquantile(med,.975)})
 curve=pd.DataFrame(out);curve.to_csv(O/'R1_SEM_sensitivity_curve.csv',index=False);critical=float(curve.loc[curve.any_Holm_p_lt_0_05,'SEM'].max()) if curve.any_Holm_p_lt_0_05.any() else None;critical_raw=float(curve.loc[curve.any_unadjusted_p_lt_0_05,'SEM'].max()) if curve.any_unadjusted_p_lt_0_05.any() else None
 cd8=x[x.celltype.eq('CD8T')][['metric','stratum','n_donors','SEM','SEM_CI_low','SEM_CI_high']].copy();cd8['critical_SEM_Holm']=critical
 for c in ['SEM','SEM_CI_low','SEM_CI_high']:cd8['critical_over_'+c]=critical/cd8[c] if critical else np.nan;cd8[c+'_over_critical']=cd8[c]/critical if critical else np.nan
 cd8.to_csv(O/'R1_critical_SEM_vs_CD8_estimates.csv',index=False);summary={'status':'PASS','seed':SEED,'n_grid_points':len(curve),'critical_SEM_Holm_grid':critical,'critical_SEM_unadjusted_grid_sensitivity':critical_raw,'critical_definition':'Largest grid point with either two-sided empirical statistic Holm(2)-adjusted p<0.05; no interpolation','n_simulations_per_grid_point':B,'common_random_numbers':'Same100000x38 standard normals across grid; pointwise tests, not a multiplicity-adjusted confidence band for all56 grid points','null_center':'observed full-cohort mean, not zero','not_a_test_of_zero_average_BCG_effect':True,'negative_median_empty_simulations':'Undefined when zero negative changes; exclude only from median p as in round2, record valid count; count p uses all simulations','all_residual_CD8_SEM_lower_CIs_above_critical':bool(cd8[cd8.metric.eq('age_acceleration')].SEM_CI_low.gt(critical).all()) if critical else None,'recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(O/'R1_summary.json').write_text(json.dumps(summary,indent=2));print('R1',summary,flush=True)
if __name__=='__main__':main()
