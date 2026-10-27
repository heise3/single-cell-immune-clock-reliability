#!/usr/bin/env python3
"""Prespecified CD8-only RTM null; 90d primary and 7-9d sensitivity."""
from pathlib import Path
import sys,json,hashlib,zipfile,io,re
import numpy as np,pandas as pd
from scipy import stats
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';SEED=20261025;B=100000
sys.path.insert(0,str(O));from run_G1_prime import rng_for

def ep(obs,null):
 null=np.asarray(null);null=null[np.isfinite(null)];return min(1.,2*min((np.sum(null<=obs)+1)/(len(null)+1),(np.sum(null>=obs)+1)/(len(null)+1)))
def holm(p):
 p=np.array(p);order=np.argsort(p);z=np.maximum.accumulate(p[order]*np.arange(len(p),0,-1));out=np.empty(len(p));out[order]=np.minimum(z,1);return out

def main():
 gate=pd.read_csv(O/'G1_prime_results.csv');cd8=gate[gate.celltype.eq('CD8T')];assert len(cd8)==1
 if cd8.iloc[0].G1_prime_status!='PASS':
  (O/'RTM_status.json').write_text(json.dumps({'status':'NOT_EVALUABLE','reason':'CD8 failed G1 prime; no substitute analysis'},indent=2));return
 reliability=pd.read_csv(O/'reliability_results.csv');src=P/'work/scImmuAging_source_data/SD_F5/CD8T.txt';d=pd.read_csv(src,sep='\t');d=d[d.condition.isin(['T0_RPMI','T3m_RPMI'])];g=d.groupby(['donor_id','condition']).agg(prediction=('Prediction','mean'),age=('age','first'),n_pseudocells=('Prediction','size')).reset_index();p=g.pivot(index='donor_id',columns='condition',values='prediction');p.columns=['baseline','followup'] if list(p.columns)==['T0_RPMI','T3m_RPMI'] else p.columns
 assert len(p)==38 and p.notna().all().all();p['change']=p.followup-p.baseline;p['age']=g[g.condition.eq('T0_RPMI')].set_index('donor_id').age;p.to_csv(O/'BCG_38_paired_CD8.csv')
 delta=p.change.to_numpy();mu=delta.mean();observed_n=int((delta<0).sum());observed_med=float(np.median(delta[delta<0]));centered=delta-np.median(delta)
 # Symmetrized sign-flip skewness diagnostic. This is approximate because
 # the center is estimated by the sample median, not specified a priori.
 observed_skew=float(stats.skew(centered,bias=False));rng=rng_for('BCG_symmetry');sk=[]
 for _ in range(100):
  a=centered[None,:]*rng.choice([-1,1],size=(1000,len(delta)));sk.extend(stats.skew(a,axis=1,bias=False))
 symp=(np.sum(np.abs(sk)>=abs(observed_skew))+1)/(B+1)
 summary={'n_paired_donors':38,'mean_change':mu,'median_change':float(np.median(delta)),'n_negative':observed_n,'negative_subgroup_median':observed_med,'symmetry_skewness':observed_skew,'signflip_symmetry_p':symp,'symmetry_caveat':'Exploratory sign-flip skewness diagnostic around sample median; failure to reject does not establish symmetry.','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source_aggregation':'mean of native 100 pseudocells per donor/condition, RPMI only; A1 median applies only to G1 prime','causal_boundary':'Paired changes without untreated controls do not identify a causal BCG effect.'};rows=[];corr=[]
 for stratum,role in [('75-105d','PRIMARY'),('7-9d','SENSITIVITY')]:
  rr=reliability[(reliability.clock=='sc-ImmuAging')&(reliability.celltype=='CD8T')&(reliability.metric=='age_acceleration')&(reliability.stratum==stratum)];assert len(rr)==1;sem=float(rr.iloc[0].SEM);a=rng_for('RTM',stratum).normal(mu,np.sqrt(2)*sem,size=(B,len(delta)));n=(a<0).sum(axis=1);neg=np.where(a<0,a,np.nan);med=np.nanmedian(neg,axis=1);pn=ep(observed_n,n);pm=ep(observed_med,med);adj=holm([pn,pm]);np.savez_compressed(O/('RTM_null_distributions_'+stratum+'.npz'),negative_count=n,negative_subgroup_median=med);rows.append({'stratum':stratum,'role':role,'SEM_CD8_age_acceleration':sem,'null_mean':mu,'null_SD':np.sqrt(2)*sem,'observed_n_negative':observed_n,'null_n_median':np.median(n),'null_n_p025':np.quantile(n,.025),'null_n_p975':np.quantile(n,.975),'count_empirical_p':pn,'count_Holm_p':adj[0],'observed_negative_median':observed_med,'null_negative_median_p025':np.nanquantile(med,.025),'null_negative_median_p975':np.nanquantile(med,.975),'negative_median_empirical_p':pm,'negative_median_Holm_p':adj[1],'null_rejected':bool((adj<.05).any()),'n_simulations':B,'n_valid_negative_medians':int(np.isfinite(med).sum())})
  x=p.baseline.to_numpy();v=np.var(x,ddof=1);cov=np.cov(x,delta,ddof=1)[0,1];den=v-sem**2;slope=(cov+sem**2)/den if den>0 else np.nan
  corr.append({'stratum':stratum,'baseline_change_Pearson_r':stats.pearsonr(x,delta).statistic,'baseline_change_Pearson_p':stats.pearsonr(x,delta).pvalue,'observed_slope':cov/v,'SEM2':sem**2,'Blomqvist_corrected_baseline_variance':den,'Blomqvist_corrected_slope':slope,'Blomqvist_status':'ESTIMABLE_UNDER_CLASSICAL_INDEPENDENT_ERROR_ASSUMPTIONS' if den>0 else 'NOT_IDENTIFIABLE_NONPOSITIVE_TRUE_BASELINE_VARIANCE'})
 import statsmodels.api as sm
 fit=sm.OLS(p.followup,sm.add_constant(p[['baseline','age']])).fit();fit.summary().as_text();(O/'BCG_ANCOVA.txt').write_text(fit.summary().as_text()+'\nNo untreated arm: this regression cannot estimate a causal intervention effect.\n')
 summary['baseline_change_correlation_is_mathematically_coupled']=True
 age_wide=g.pivot(index='donor_id',columns='condition',values='age');summary['source_age_equal_before_after']=bool(age_wide['T0_RPMI'].eq(age_wide['T3m_RPMI']).all());summary['delta_scale_note']='With identical released age before and after, any common intercept+age adjustment cancels in the paired difference. External K0 SEM transportability remains an assumption.'
 inventory=[];candidates=[];source_zip=P/'work/sources/43587_2025_819_MOESM6_ESM.zip'
 with zipfile.ZipFile(source_zip) as zf:
  for name in zf.namelist():
   if not name.lower().endswith(('.txt','.tsv','.csv')) or '__MACOSX' in name:continue
   try:
    tab=pd.read_csv(io.BytesIO(zf.read(name)),sep=None,engine='python');ifn=[col for col in tab if re.search(r'ifn|interferon',str(col),re.I)];don=[col for col in tab if re.fullmatch(r'donor[_ .]?id|donor|subject[_ .]?id',str(col),re.I)];inventory.append({'file':name,'columns':'|'.join(map(str,tab.columns)),'n_rows':len(tab),'IFN_columns':'|'.join(map(str,ifn)),'donor_columns':'|'.join(map(str,don))})
    if ifn and don:candidates.append((name,tab,don[0],ifn))
   except (UnicodeDecodeError,pd.errors.ParserError):inventory.append({'file':name,'columns':'UNPARSEABLE_TEXT_TABLE'})
 pd.DataFrame(inventory).to_csv(O/'Figure5_source_table_inventory.csv',index=False)
 summary['IFN_status']='NOT_EVALUABLE_SOURCE_HAS_NO_DONOR_MATCHED_IFN_SCORE';summary['IFN_missing']='Individual baseline IFN score matched to the same38 donor IDs. No score is inferred from condition labels or group summaries.'
 if candidates:
  ifnrows=[]
  for name,tab,don,columns in candidates:
   cond=next((col for col in tab if str(col).lower() in ['condition','visit','timepoint']),None)
   if cond is not None:tab=tab[tab[cond].astype(str).str.contains(r'T0|baseline',case=False,regex=True)]
   for col in columns:
    base=tab.assign(IFN_numeric=pd.to_numeric(tab[col],errors='coerce')).groupby(don).IFN_numeric.median();joined=p.join(base.rename('baseline_IFN'),how='inner').dropna(subset=['baseline_IFN']);fit_ifn=sm.OLS(joined.followup,sm.add_constant(joined[['baseline','age','baseline_IFN']])).fit(cov_type='HC3') if len(joined)>=10 and joined.baseline_IFN.nunique()>1 else None
    if fit_ifn is not None:ifnrows.append({'source_file':name,'score_column':col,'n_donors':len(joined),'adjusted_followup_IFN_beta':fit_ifn.params['baseline_IFN'],'HC3_p':fit_ifn.pvalues['baseline_IFN'],'HC3_CI_low':fit_ifn.conf_int().loc['baseline_IFN',0],'HC3_CI_high':fit_ifn.conf_int().loc['baseline_IFN',1],'model':'followup clock ~ baseline clock + chronological age + baseline IFN; donor medians for repeated score rows'})
  if ifnrows:pd.DataFrame(ifnrows).to_csv(O/'BCG_IFN_baseline_adjusted.csv',index=False);summary['IFN_status']='ANALYZED_DONOR_MATCHED_PUBLISHED_SCORE'
  else:summary['IFN_status']='NOT_EVALUABLE_IFN_COLUMNS_NOT_MATCHABLE_TO_PAIRED_BASELINE_DONORS'

 pd.DataFrame(rows).to_csv(O/'RTM_null_results.csv',index=False);pd.DataFrame(corr).to_csv(O/'RTM_Blomqvist.csv',index=False);(O/'RTM_summary.json').write_text(json.dumps(summary,indent=2));print(summary);print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
