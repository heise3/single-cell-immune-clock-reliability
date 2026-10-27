#!/usr/bin/env python3
"""Polite documented bulk-search endpoint with Boolean query refinement."""
from pathlib import Path
import json,time,requests,hashlib,datetime,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';D=P/'work/round3/novelty_receipts'
QUERIES=['transcriptomic +clock +(reliability | "test-retest")','"single-cell" +aging +clock +ICC','TranslAGE +transcriptomic','"immune aging" +clock +(responders | "regression to the mean")']
def main():
 logs=[];records=[];url='https://api.semanticscholar.org/graph/v1/paper/search/bulk'
 for i,q in enumerate(QUERIES):
  params={'query':q,'year':'2025-2026','fields':'title,year,publicationDate,authors,externalIds,url,abstract'};row={'engine':'SemanticScholar_bulk','query':q,'url':url,'parameters':params,'queried_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};attempt=[];data={}
  for k in range(2):
   try:
    r=requests.get(url,params=params,timeout=45);f=D/f'SemanticScholar_bulk_{i+1}_{k+1}.json';f.write_bytes(r.content);attempt.append({'http_status':r.status_code,'receipt':str(f),'sha256':hashlib.sha256(r.content).hexdigest()})
    if r.status_code==200:data=r.json();break
    if r.status_code==429 and k==0:time.sleep(min(30,max(5,int(r.headers.get('Retry-After','30')))));continue
    break
   except (requests.RequestException,ValueError) as e:attempt.append({'error':str(e)})
  row.update(attempts=attempt,status='RESULTS_RETRIEVED' if 'data' in data else 'UNKNOWN_EXTERNAL_API_ACCESS',total_hits=data.get('total'),returned_records=len(data.get('data',[])),bounded_one_bulk_page=True);logs.append(row)
  for h in data.get('data',[]):records.append({'engine':'SemanticScholar_bulk','query':q,'id':h['paperId'],'title':h.get('title'),'year':h.get('year'),'date':h.get('publicationDate'),'doi':(h.get('externalIds') or {}).get('DOI'),'url':h.get('url'),'abstract':h.get('abstract') or ''})
  time.sleep(2)
 (O/'R3_SemanticScholar_bulk_receipts.json').write_text(json.dumps(logs,indent=2));t=pd.DataFrame(records)
 if len(t):
  t.to_csv(D/'SemanticScholar_bulk_abstracts_INTERNAL.csv',index=False);t['abstract_sha256']=t.abstract.map(lambda x:hashlib.sha256(x.encode()).hexdigest());t['abstract_characters']=t.abstract.str.len();t.drop(columns='abstract').to_csv(O/'R3_SemanticScholar_bulk_records.csv',index=False)
 print(json.dumps([{k:v for k,v in r.items() if k not in ['parameters','attempts']} for r in logs],indent=2))
if __name__=='__main__':main()
