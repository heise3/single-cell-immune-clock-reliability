#!/usr/bin/env python3
"""K3: deterministic structural audit of public GEO metadata, no expression analysis.

Eligible visits: the explicitly non-vaccination Immune Variation series only.
Known COVID vaccine dates falling in [first,last] invalidate that pair.
No unstated biological washout period is invented. Missing COVID vaccination
dates remain unknown and are reported separately from known clean intervals.
Seed 20261025 is reserved; this metadata enumeration uses no random sampling.
"""
import csv,gzip,itertools,json,re,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parent
text=gzip.open(ROOT.parent/'work/sources/GSE271896_family.soft.gz','rt').read()
samples=[]
for block in text.split('^SAMPLE = ')[1:]:
 r={'gsm':block.splitlines()[0].strip()}
 for line in block.splitlines():
  if line.startswith('!Sample_characteristics_ch1 = '):
   key,value=line.split(' = ',1)[1].split(': ',1);r[key]=value
  if line.startswith('!Sample_supplementary_file_1 = '):r['count_url']=line.split(' = ',1)[1].replace('ftp://','https://')
 samples.append(r)
fields=sorted(set().union(*(r.keys() for r in samples)))
with (ROOT/'sample_metadata.csv').open('w') as f:
 writer=csv.DictWriter(f,fields);writer.writeheader();writer.writerows(samples)
donors=collections.defaultdict(list)
for r in samples:donors[r['subject id']].append(r)
pairs=[];summary=[]
for donor,rows in sorted(donors.items()):
 base=[r for r in rows if r['visit'].startswith('Immune Variation') and r.get('sample bridging_control')=='FALSE']
 vaccines=sorted(set(int(r[k]) for r in rows for k in r if k.startswith('covid vax_dose_') and r[k] not in ['NA','','N/A']))
 flu_days=sorted(set(int(r['days since_first_visit']) for r in rows if re.fullmatch(r'Flu Year [0-9]+ Day 0',r['visit'])))
 for a,b in itertools.combinations(sorted(base,key=lambda r:int(r['days since_first_visit'])),2):
  lo,hi=int(a['days since_first_visit']),int(b['days since_first_visit'])
  events=[v for v in vaccines if lo<=v<=hi]
  pairs.append({'donor':donor,'gsm1':a['gsm'],'gsm2':b['gsm'],'visit1':a['visit'],'visit2':b['visit'],'day1':lo,'day2':hi,'interval_days':hi-lo,'covid_vaccine_days':json.dumps(vaccines),'vaccine_in_interval':bool(events),'vaccine_timeline_status':'PARTIAL_RECORDED' if vaccines else 'NO_DATE_RECORDED','flu_vaccine_days':json.dumps(flu_days),'any_recorded_vaccine_in_interval':any(lo<=v<=hi for v in vaccines+flu_days)})
 summary.append({'donor':donor,'n_immune_variation_visits':len(base),'n_pairs':len(base)*(len(base)-1)//2,'covid_vaccine_days':vaccines})
if pairs:
 with (ROOT/'candidate_visit_pairs.csv').open('w') as f:
  w=csv.DictWriter(f,list(pairs[0]));w.writeheader();w.writerows(pairs)
out={'seed':20261025,'samples':len(samples),'donors':len(donors),'visit_counts':dict(collections.Counter(r['visit'] for r in samples)),
 'donors_ge2_immune_variation_visits':sum(r['n_immune_variation_visits']>=2 for r in summary),
 'donors_ge3_immune_variation_visits':sum(r['n_immune_variation_visits']>=3 for r in summary),
 'donors_with_any_pair_without_recorded_vaccine_in_interval':len({r['donor'] for r in pairs if not r['vaccine_in_interval']}),
 'donors_with_any_pair_without_recorded_vaccine_in_interval_and_at_least_one_vaccine_date':len({r['donor'] for r in pairs if not r['vaccine_in_interval'] and r['vaccine_timeline_status']=='PARTIAL_RECORDED'}),
 'donors_with_no_covid_vaccine_date':sum(not r['covid_vaccine_days'] for r in summary),
 'donors_with_any_pair_without_recorded_COVID_or_flu_vaccine_in_interval':len({r['donor'] for r in pairs if not r['any_recorded_vaccine_in_interval']}),
 'pairs_by_interval_days':dict(collections.Counter(r['interval_days'] for r in pairs)),
 'caveat':'No recorded event is not proof of no event. No biological washout period has been imposed. Subject status and recorded visit perturbations must be checked before K4.',
 'donor_summary':summary}
(ROOT/'K3_structure.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
selected_path=ROOT/'K4_selected_pairs.csv'
if selected_path.exists():
 selected=list(csv.DictReader(selected_path.open()));inside=[];spacing=[]
 for row in selected:
  flu=sorted({int(r['days since_first_visit']) for r in donors[row['donor']] if re.fullmatch(r'Flu Year [0-9]+ Day 0',r['visit'])})
  lo,hi=int(row['day1']),int(row['day2'])
  if any(lo<=v<=hi for v in flu):inside.append(row['donor'])
  earlier=[v for v in flu if v<lo]
  if earlier:spacing.append(lo-max(earlier))
 qa={'donors_with_any_pair_without_recorded_COVID_or_flu_vaccine_in_interval':out['donors_with_any_pair_without_recorded_COVID_or_flu_vaccine_in_interval'],'selected_pairs_with_known_flu_event_in_interval':inside,'selected_pair_days_since_last_documented_flu_vaccine_range':[min(spacing),max(spacing)]}
 (ROOT/'K3_all_recorded_vaccine_QA.json').write_text(json.dumps(qa,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='donor_summary'},ensure_ascii=False,indent=2))
