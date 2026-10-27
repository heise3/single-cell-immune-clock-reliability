#!/usr/bin/env python3
"""Recreate K0 CMV/sex extraction from the archived GEO SOFT, without inference."""
from pathlib import Path
import gzip,json
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5'
def main():
 src=P/'work/sources/GSE271896_family.soft.gz';records=[];r=None;keys=set()
 with gzip.open(src,'rt') as h:
  for ln in h:
   if ln.startswith('^SAMPLE = '):
    if r is not None:records.append(r)
    r={'gsm':ln.rstrip().split(' = ',1)[1]}
   elif r is not None and ln.startswith('!Sample_characteristics_ch1 = '):
    z=ln.rstrip().split(' = ',1)[1]
    if ': ' in z:
     k,v=z.split(': ',1);r[k]=v;keys.add(k)
 if r is not None:records.append(r)
 allm=pd.DataFrame(records);k=pd.read_csv(P/'outputs/round2_A1/provenance_inputs/K0_primary_with_released_age.csv');v=k[['gsm','donor']].merge(allm,on='gsm',validate='one_to_one');assert len(v)==320 and v.donor.nunique()==68 and v.donor.eq(v['subject id']).all();assert v.groupby('donor')[['Sex','cmv status']].nunique().max().eq(1).all();v.to_csv(O/'GEO_K0_covariate_visits.csv',index=False);d=v.groupby('donor',sort=True).agg(Sex=('Sex','first'),**{'cmv status':('cmv status','first')},n_visits=('gsm','size')).reset_index();d.to_csv(O/'GEO_K0_covariate_donors.csv',index=False);assert not any(any(s in x.lower() for s in ['bmi','weight','height','body mass']) for x in keys);(O/'GEO_covariate_schema.json').write_text(json.dumps({'all_characteristic_keys':sorted(keys),'BMI_available':False,'source':str(src.relative_to(P)),'n_K0_donors':68,'n_K0_visits':320},indent=2));print(d.groupby(['Sex','cmv status']).size().to_string())
if __name__=='__main__':main()
