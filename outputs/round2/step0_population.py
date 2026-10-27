#!/usr/bin/env python3
"""Freeze the user-designated HIHA K0 population and compare round 1."""
from pathlib import Path
import hashlib,json,shutil,itertools
import pandas as pd
ROOT=Path(__file__).resolve().parent;PROJECT=ROOT.parents[1];WORK=PROJECT/'work/round2'
import os
SOURCE=Path(os.environ.get('HIHA_K0_SOURCE', str(PROJECT/'cohort_audit/outputs/HIHA_gate_audit')))
SNAP=WORK/'K0_snapshot';SNAP.mkdir(exist_ok=True,parents=True)
files=['K0_v2/sample_eligibility.csv','K0_v2/donor_eligibility.csv','K0_v2/window_summary.json','K0_donor_visits.csv','K0_sample_audit.csv','K0_metadata_summary.json','K0_CELLxGENE_assets.json','GEO_sample_metadata.csv','gate_contract.json']
receipts=[]
for name in files:
 p=SOURCE/name;dst=SNAP/name;dst.parent.mkdir(exist_ok=True,parents=True);shutil.copyfile(p,dst)
 receipts.append({'source':str(p),'snapshot':str(dst.relative_to(PROJECT)),'sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'bytes':dst.stat().st_size})
(ROOT/'K0_source_manifest.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2))
d=pd.read_csv(SNAP/'K0_v2/donor_eligibility.csv',keep_default_na=False)
s=pd.read_csv(SNAP/'K0_v2/sample_eligibility.csv',keep_default_na=False)
audit=pd.read_csv(SNAP/'K0_sample_audit.csv',keep_default_na=False)
main_d=d[(d.window_days==90)&d.donor_ge4.astype(str).str.lower().eq('true')]
main_ids=set(main_d.donor)
main=s[(s.window_days==90)&s.eligible_visit.astype(str).str.lower().eq('true')&s.donor.isin(main_ids)].copy()
assert len(main_ids)==68 and len(main)==320
main=main.merge(audit[['gsm','batch id','birth year','Sex','cmv status','specimen.specimenGuid','dataset_id','dataset_version_id']],on='gsm',validate='one_to_one')
main.to_csv(ROOT/'K0_primary_68_donors_320_visits.csv',index=False)
pd.DataFrame(sorted(main_ids),columns=['donor']).to_csv(ROOT/'K0_primary_donors.csv',index=False)
batch=main.groupby(['batch id','visit']).agg(n_samples=('gsm','nunique'),n_donors=('donor','nunique')).reset_index()
batch.to_csv(ROOT/'K0_primary_batch_ledger.csv',index=False)
pairs=pd.read_csv(PROJECT/'outputs/candidate_visit_pairs.csv',keep_default_na=False)
r1_89=set(pairs.donor)
r1_87=set(pairs.loc[pairs.any_recorded_vaccine_in_interval.astype(str).str.lower().eq('false'),'donor'])
assert len(r1_89)==89 and len(r1_87)==87
all_donors=sorted(set(d.donor));lookup=d[d.window_days==90].set_index('donor');rows=[]
for donor in all_donors:
 row=lookup.loc[donor];reason=[]
 if not bool(row.complete_covid_dates):reason.append('K0 excludes incomplete COVID dose 1/2 dates')
 if bool(row.ambiguous_visit_labels):reason.append('K0 excludes ambiguous visit labels')
 if row.n_distinct_eligible_days<4:reason.append(f'K0 has {row.n_distinct_eligible_days} eligible visits (<4) after frozen 90-day exclusions')
 if donor not in r1_89:reason.append('Round1 fewer than two Immune Variation visits')
 elif donor not in r1_87:reason.append('Round1 all candidate IV pairs contain recorded vaccination')
 if donor in main_ids and donor in r1_87:reason.append('Included in both; visit sets and washout definitions still differ')
 rows.append({'donor':donor,'K0_primary_68':donor in main_ids,'round1_89':donor in r1_89,'round1_87':donor in r1_87,'K0_eligible_visit_count':int(row.n_distinct_eligible_days),'complete_COVID_dates':row.complete_covid_dates,'ambiguous_visit_labels':row.ambiguous_visit_labels,'reason':'; '.join(reason)})
comparison=pd.DataFrame(rows);comparison.to_csv(ROOT/'population_comparison_all96.csv',index=False)
comparison[(comparison.K0_primary_68!=comparison.round1_89)|(comparison.K0_primary_68!=comparison.round1_87)].to_csv(ROOT/'population_difference_donors.csv',index=False)
# Pre-specified nearest single pair per donor and stratum.
strata=[('7-9d',7,9,7),('75-105d',75,105,90),('330-400d',330,400,365)]
pr=[]
for donor,grp in main.groupby('donor',sort=True):
 candidates=[]
 for (_,a),(_,b) in itertools.combinations(grp.sort_values(['relative_day','gsm']).iterrows(),2):
  delta=int(b.relative_day)-int(a.relative_day)
  for label,lo,hi,target in strata:
   if lo<=delta<=hi:candidates.append({'donor':donor,'stratum':label,'gsm1':a.gsm,'gsm2':b.gsm,'day1':int(a.relative_day),'day2':int(b.relative_day),'interval_days':delta,'target_distance':abs(delta-target)})
 for label,lo,hi,target in strata:
  opts=[p for p in candidates if p['stratum']==label]
  if opts:pr.append(min(opts,key=lambda p:(p['target_distance'],p['day1'],p['gsm1'],p['gsm2'])))
pd.DataFrame(pr).to_csv(ROOT/'K0_primary_pair_ledger.csv',index=False)
summary={'seed':20261025,'primary_donors':68,'primary_visits':320,'primary_batches':int(main['batch id'].nunique()),'round1_89_overlap':len(main_ids&r1_89),'round1_87_overlap':len(main_ids&r1_87),'K0_only_vs87':sorted(main_ids-r1_87),'round1_87_only':sorted(r1_87-main_ids),'pair_donors_by_stratum':pd.Series([p['stratum'] for p in pr]).value_counts().to_dict(),'main_window':'Frozen HIHA K0 90-day ledger; no eligibility redefinition','sensitivity':'Round1 87 donor cohort only','status':'PASS'}
(ROOT/'step0_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2))
