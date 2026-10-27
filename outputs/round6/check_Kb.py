#!/usr/bin/env python3
"""Biological sampling units and nominal annual-wave intervals; no clock values."""
from pathlib import Path
import pandas as pd,numpy as np,json,datetime,re
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6'
def wave(tube):
 # F prefixes female IDs (e.g.FA02->FG02),not an extra temporal wave.
 s=re.sub(r'^F(?=[A-Z])','',tube);c=re.match(r'^([A-Z])\d+$',s);assert c,tube
 return {**dict.fromkeys('ABCDE',0),**dict.fromkeys('GHJKL',1),**dict.fromkeys('MNPST',2)}[c[1]]
def main():
 assert json.loads((O/'K_a_status.json').read_text())['status']=='PASS';z=pd.read_csv(P/'work/round6/extracted/b_cells_metadata.csv',index_col=0);v=pd.read_csv(O/'K_b_B_donor_tube_metadata_screen.csv');assert v.n_age.eq(1).all() and not v[['Donor_id','Tube_id']].isna().any().any();v['visit_wave']=v.Tube_id.map(wave);assert not v.duplicated(['Donor_id','visit_wave']).any();rows=[];allpairs=[]
 for donor,g in v.groupby('Donor_id'):
  g=g.sort_values('visit_wave');assert g.Age_min.is_monotonic_increasing
  if len(g)<2:continue
  a,b=g.iloc[0],g.iloc[-1];row={'donor':donor,'visit1':a.Tube_id,'visit2':b.Tube_id,'wave1':a.visit_wave,'wave2':b.visit_wave,'wave_gap':b.visit_wave-a.visit_wave,'interval_nominal_years':b.visit_wave-a.visit_wave,'interval_actual_days':np.nan,'integer_age_difference_NOT_actual_interval':b.Age_min-a.Age_min,'age1':a.Age_min,'age2':b.Age_min,'n_available_visits':len(g),'same_source_B_batch':a.Batches==b.Batches,'batch1_B':a.Batches,'batch2_B':b.Batches};rows.append(row)
  for i in range(len(g)):
   for j in range(i+1,len(g)):allpairs.append({'donor':donor,'visit1':g.iloc[i].Tube_id,'visit2':g.iloc[j].Tube_id,'wave_gap':g.iloc[j].visit_wave-g.iloc[i].visit_wave,'actual_days':np.nan})
 pair=pd.DataFrame(rows);pair.to_csv(O/'K_b_preQC_earliest_latest_pairs.csv',index=False);pd.DataFrame(allpairs).to_csv(O/'K_b_all_within_donor_pairs.csv',index=False);v.to_csv(O/'K_b_donor_visit_wave_ledger.csv',index=False);dist=pair.groupby('interval_nominal_years').size().reset_index(name='n_donors');dist['basis']='nominalannual-wavegap;not exactdatesor integerage subtraction';dist['actual_day_interval_available']=False;dist.to_csv(O/'K_b_nominal_interval_distribution.csv',index=False)
 # Same-source time-design evidence: officialdataset3-yearcollection plus an
 # independently published reanalysis explicitly describing annual collection.
 status={'gate':'K-b','status':'PASS','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_source_donors':166,'n_source_biological_visit_units':317,'n_donors_ge2visits':len(pair),'visits_per_donor_distribution':v.groupby('Donor_id').size().value_counts().sort_index().to_dict(),'nominal_earliest_latest_year_interval_distribution':dict(zip(dist.interval_nominal_years.astype(str),dist.n_donors.astype(int))),'actual_day_interval_distribution':'NOT_RELEASED;no dates or time offsets in sourceCSV. Integerage differences retained solely asmetadataQC,not used to derive intervals.','wave_mapping':'baselineABCDE,followupGHJKL,secondfollowupMNPST;optionalF female-prefix stripped;all visitwave/donor pairsunique,age monotone. Wave decoding is an explicitly labelled inference from source identifiers consistent with annual design,not a released exact date.','health':'HCA specimen and donor disease=normal;healthy adult study25-85years,unperturbed longitudinal samples;not vaccination/recovery visits.','time_design_source':'OriginalHCA/AWSdescribe3-yearcollection;explicitannualcollection in primaryreanalysis Methods:Terekhova2023dataset section ofResearchSquare rs-6784948,subsequently Sex-specific trajectories of nonlinear immune aging at single-cell level;exactdatesnotavailable.','interval_interpretation':'Pooled reliability over nominal1/2-year earliest-latest visits only;not exact90day/365day strata.','technical_library_collapse':'Tube_id is visit;File_name technical libraries within Tube will be pooled;305Bvisits2libraries,6one,6three;never count libraries as visits.'};assert len(pair)>=15;(O/'K_b_status.json').write_text(json.dumps(status,indent=2));f=P/'outputs/gate_ledger.json';q=json.loads(f.read_text());q['round6']['K_b_result']=status;f.write_text(json.dumps(q,ensure_ascii=False,indent=2));print(json.dumps(status,indent=2))
if __name__=='__main__':main()
