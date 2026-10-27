#!/usr/bin/env python3
"""Anonymous API queries; preserve status and bounds, never treat429 as no hit."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests,json,time,hashlib,datetime
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';D=P/'work/round3/novelty_receipts';D.mkdir(exist_ok=True);Q=['transcriptomic clock reliability','transcriptomic clock test-retest','single-cell aging clock ICC','TranslAGE transcriptomic','immune aging clock responders regression to the mean']
def run(engine):
 logs=[];records=[]
 for i,q in enumerate(Q):
  url='https://api.openalex.org/works' if engine=='OpenAlex' else 'https://api.semanticscholar.org/graph/v1/paper/search'
  params={'search':q,'filter':'from_publication_date:2025-01-01,to_publication_date:2026-10-26','per-page':100,'select':'id,doi,title,publication_year,publication_date,authorships,abstract_inverted_index,primary_location,open_access'} if engine=='OpenAlex' else {'query':q,'year':'2025-2026','limit':100,'fields':'title,year,publicationDate,authors,externalIds,url,abstract,openAccessPdf'}
  row={'engine':engine,'query':q,'url':url,'parameters':params,'queried_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};attempts=[];data={}
  for attempt in range(3):
   try:
    r=requests.get(url,params=params,headers={'User-Agent':'single-cell-clock-reliability-audit/3'},timeout=45);f=D/f'{engine}_{i+1}_attempt{attempt+1}.json';f.write_bytes(r.content);attempts.append({'http_status':r.status_code,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest(),'receipt':str(f)});row['request_url']=r.url
    if r.status_code==200:data=r.json();break
    if r.status_code not in [429,500,502,503,504]:break
    if attempt<2:time.sleep([5,30][attempt])
   except (requests.RequestException,ValueError) as e:attempts.append({'error':str(e)})
  hits=data.get('results',[]) if engine=='OpenAlex' else data.get('data',[]);row.update(attempts=attempts,status='RESULTS_RETRIEVED' if data and ('results' in data or 'data' in data) else 'UNKNOWN_EXTERNAL_API_ACCESS',total_hits=data.get('meta',{}).get('count') if engine=='OpenAlex' else data.get('total'),returned_records=len(hits),bounded_first_page=True);logs.append(row)
  for h in hits:
   if engine=='OpenAlex':
    abstract=h.get('abstract_inverted_index') or {};abstract_text=' '.join(w for _,w in sorted((pos,w) for w,pp in abstract.items() for pos in pp));rec={'engine':engine,'query':q,'id':h['id'],'title':h.get('title'),'year':h.get('publication_year'),'date':h.get('publication_date'),'doi':h.get('doi'),'url':(h.get('primary_location') or {}).get('landing_page_url'),'authors':'|'.join((z.get('author') or {}).get('display_name','') for z in h.get('authorships',[])),'abstract':abstract_text}
   else:rec={'engine':engine,'query':q,'id':h['paperId'],'title':h.get('title'),'year':h.get('year'),'date':h.get('publicationDate'),'doi':(h.get('externalIds') or {}).get('DOI'),'url':h.get('url'),'authors':'|'.join(z.get('name','') for z in h.get('authors',[])),'abstract':h.get('abstract') or ''}
   records.append(rec)
  time.sleep(2)
 return logs,records
def main():
 logs=[];records=[]
 with ThreadPoolExecutor(max_workers=2) as ex:
  for l,r in ex.map(run,['OpenAlex','SemanticScholar']):logs+=l;records+=r
 (O/'R3_API_search_receipts.json').write_text(json.dumps(logs,indent=2));pd.DataFrame([{k:(json.dumps(v) if isinstance(v,(dict,list)) else v) for k,v in r.items()} for r in logs]).to_csv(O/'R3_API_search_log.csv',index=False)
 if records:
  t=pd.DataFrame(records);t.to_csv(D/'all_records_with_abstracts_INTERNAL.csv',index=False);t['abstract_sha256']=t.abstract.map(lambda s:hashlib.sha256(s.encode()).hexdigest());t['abstract_characters']=t.abstract.str.len();t.drop(columns='abstract').to_csv(O/'R3_API_bibliographic_records.csv',index=False)
 print(json.dumps([{'engine':r['engine'],'query':r['query'],'status':r['status'],'hits':r['total_hits'],'returned':r['returned_records']} for r in logs],indent=2))
if __name__=='__main__':main()
