#!/usr/bin/env python3
"""Descriptive effect magnitudes reconstructed from published source data.

Source Fig 5: unstimulated RPMI, paired donor mean of 100 pseudocells.
Source Fig 4: the exact-age 63 subgroup named in the paper Methods; a common
age-only bias correction cancels in contrasts at this same chronological age.
These are source-data reconstructions, not numbers quoted from the article,
not a single overall BCG rejuvenation effect, and not new clinical inference.
"""
from pathlib import Path
import json,subprocess,zipfile
import pandas as pd
ROOT=Path(__file__).resolve().parent;WORK=ROOT.parent/'work';SOURCES=WORK/'sources'
base='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43587-025-00819-z/MediaObjects/'
for i in [5,6]:
 p=SOURCES/f'43587_2025_819_MOESM{i}_ESM.zip'
 if not p.exists():subprocess.run(['curl','-fsSL','--max-time','40','-o',str(p),base+p.name],check=True)
 with zipfile.ZipFile(p) as z:z.extractall(WORK/'scImmuAging_source_data')
b=pd.read_csv(WORK/'scImmuAging_source_data/SD_F5/CD8T.txt',sep='\t')
a=b.groupby(['donor_id','condition']).Prediction.mean().unstack()
delta=(a.T3m_RPMI-a.T0_RPMI).dropna()
pd.DataFrame({'donor':delta.index,'CD8T_RPMI_T3m_minus_T0_years':delta.values}).to_csv(ROOT/'published_BCG_CD8_donor_changes.csv',index=False)
ar=delta[delta<0]
c=pd.read_csv(WORK/'scImmuAging_source_data/SD_F4/MONO_bonn_fordistribution.txt',sep='\t')
c=c.groupby(['donor_id','condition','age']).Prediction.mean().reset_index()
c=c[c.age==63]
pd.DataFrame(c).to_csv(ROOT/'published_COVID_MONO_age63.csv',index=False)
means=c.groupby('condition').Prediction.mean();counts=c.groupby('condition').size()
out={'source_DOI':'10.1038/s43587-025-00819-z','BCG':{'source_url':base+'43587_2025_819_MOESM6_ESM.zip','condition':'RPMI unstimulated, 3 months minus baseline, donor means of 100 pseudocells','n_paired':len(delta),'all_donor_mean_years':delta.mean(),'all_donor_median_years':delta.median(),'all_donor_range_years':[delta.min(),delta.max()],'AR_n':len(ar),'AR_median_years':ar.median(),'AR_range_years':[ar.min(),ar.max()],'caveat':'AR is the outcome-defined rejuvenating subgroup, not the overall intervention effect.'},'COVID_MONO':{'source_url':base+'43587_2025_819_MOESM5_ESM.zip','subgroup':'Bonn cohort, exact chronological age 63, as an age-matched comparison described in Methods','n_by_condition':counts.to_dict(),'mean_prediction_by_condition_years':means.to_dict(),'mild_minus_control_years':means['mild']-means['control'],'severe_minus_control_years':means['severe']-means['control'],'caveat':'Descriptive small age-matched subgroup reconstruction from raw source predictions. It is not the average effect for all COVID cohorts. Any common age-only calibration cancels at age 63.'}}
(ROOT/'published_effect_magnitudes.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
print(json.dumps(out,ensure_ascii=False,indent=2))
