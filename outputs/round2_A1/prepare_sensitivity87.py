#!/usr/bin/env python3
"""Freeze round1 87-donor sensitivity before prediction; never replace K0."""
from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];BASE=P/'outputs/round2_A1';O=BASE/'sensitivity87';O.mkdir(exist_ok=True);W=P/'work/round2_A1/sensitivity87';(W/'cache').mkdir(parents=True,exist_ok=True)
pairs=pd.read_csv(P/'outputs/candidate_visit_pairs.csv',keep_default_na=False);pairs=pairs[pairs.any_recorded_vaccine_in_interval.astype(str).str.lower().eq('false')].copy();assert pairs.donor.nunique()==87
rows=[]
for donor,g in pairs.groupby('donor'):
 for name,lo,hi,target in [('7-9d',7,9,7),('75-105d',75,105,90)]:
  a=g[g.interval_days.between(lo,hi)].copy();a['target_distance']=(a.interval_days-target).abs()
  if len(a):
   r=a.sort_values(['target_distance','day1','gsm1','gsm2']).iloc[0].to_dict();r['stratum']=name;rows.append(r)
pd.DataFrame(rows).to_csv(O/'sensitivity87_pair_ledger.csv',index=False)
# Include every visit participating in an eligible original pair, preserving all
# 87 original donors even when a donor has no exact prespecified interval pair.
gsms=sorted(set(pairs.gsm1)|set(pairs.gsm2));m=pd.read_csv(P/'outputs/sample_metadata.csv',keep_default_na=False);m=m[m.gsm.isin(gsms)].copy();m=m.rename(columns={'subject id':'donor','days since_first_visit':'relative_day'});m['specimen.specimenGuid']=m.count_url.str.extract(r'_(PB[0-9]+-[0-9]+)_',expand=False)
a=pd.read_csv(P/'outputs/round2/released_sample_age_metadata.csv');m=m.merge(a,on=['specimen.specimenGuid'],validate='one_to_one');assert m.donor.nunique()==87 and m.gsm.nunique()==len(gsms) and m['sample.subjectAgeAtDraw'].notna().all();m.to_csv(O/'sensitivity87_visit_ledger.csv',index=False)
for name in ['HIHA_Census_dataset_records.csv','HIHA_Census_var_0a8e3443-c3e2-4918-84db-0495657d9175.csv','RAM_mapped_gene_panels.csv','gene_mapping_HIHA_Census.csv','G1_prime_results.csv']:
 shutil.copy2(BASE/name,O/name)
shutil.copy2(BASE/'run_G1_prime.py',O/'run_G1_prime.py')
for did in pd.read_csv(BASE/'HIHA_Census_dataset_records.csv').dataset_id:
 for stem in ['HIHA_Census_obs_']:
  dst=W/'cache'/f'{stem}{did}.parquet';src=P/'work/round2_A1/cache'/dst.name
  if not dst.exists() and not dst.is_symlink():dst.symlink_to(src)
for name in ['fetch_HIHA_metadata.py','fetch_main_Census_counts.py','score_main_clocks.py','compute_reliability.py','run_primary_pipeline.py']:
 s=(BASE/name).read_text().replace('P=Path(__file__).resolve().parents[2]','P=Path(__file__).resolve().parents[3]').replace("O=P/'outputs/round2_A1'","O=P/'outputs/round2_A1/sensitivity87'").replace("W=P/'work/round2_A1'","W=P/'work/round2_A1/sensitivity87'")
 s=s.replace("W=P/'work/round2_A1/cache'","W=P/'work/round2_A1/sensitivity87/cache'")
 s=s.replace("P/'outputs/round2/K0_primary_with_released_age.csv'","O/'sensitivity87_visit_ledger.csv'").replace("P/'outputs/round2/K0_primary_pair_ledger.csv'","O/'sensitivity87_pair_ledger.csv'")
 s=s.replace('a.gsm.nunique()==320 and a.donor.nunique()==68',f'a.gsm.nunique()=={len(m)} and a.donor.nunique()==87').replace('obs.gsm.nunique()==320 and obs.donor.nunique()==68',f'obs.gsm.nunique()=={len(m)} and obs.donor.nunique()==87').replace("'n_donors':68,'n_visits':320",f"'n_donors':87,'n_visits':{len(m)}")
 # Reuse the same fully downloaded public metadata segments; these contain no X.
 s=s.replace("P/'work/round2_A1/sensitivity87/HIHA_metadata_segments'","P/'work/round2_A1/HIHA_metadata_segments'")
 (O/name).write_text(s)
report={'seed':20261025,'population':'round1 87 donor sensitivity only','n_donors':87,'n_visits':len(m),'original_eligible_pairs':len(pairs),'prespecified_interval_pairs':pd.DataFrame(rows).stratum.value_counts().to_dict(),'no_approximately_one_year_IV_pairs':True,'eligibility':'original pair has no recorded COVID or influenza vaccination within interval; original incomplete timelines retained as sensitivity, not claimed K0 steady-state','primary_remains_K0_68_320':True};(O/'sensitivity87_contract.json').write_text(json.dumps(report,indent=2));print(report)
