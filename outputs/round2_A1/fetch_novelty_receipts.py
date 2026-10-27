#!/usr/bin/env python3
"""Public primary-source and PubMed receipts. Access blocks remain explicit."""
from pathlib import Path
import requests,json,datetime,hashlib,re
from html.parser import HTMLParser
from concurrent.futures import ThreadPoolExecutor
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';C=P/'work/round2_A1/novelty_receipts';C.mkdir(exist_ok=True)
SOURCES={'IDEAL_Age':'https://link.springer.com/article/10.1186/s13059-026-04188-7','Borrus_2026':'https://onlinelibrary.wiley.com/doi/10.1111/acel.70756','Sehgal_2026':'https://onlinelibrary.wiley.com/doi/full/10.1111/acel.70635'}
class Plain(HTMLParser):
 def __init__(self):super().__init__();self.parts=[];self.links=[];self.skip=0
 def handle_starttag(self,t,a):
  if t in ['script','style']:self.skip+=1
  if t=='a':
   href=dict(a).get('href')
   if href:self.links.append(href)
 def handle_endtag(self,t):
  if t in ['script','style']:self.skip=max(0,self.skip-1)
 def handle_data(self,d):
  if not self.skip:self.parts.append(d)
def fetch(item):
 name,url=item;row={'name':name,'url':url,'queried_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  r=requests.get(url,timeout=45);file=C/(name+'.html');file.write_bytes(r.content);row.update(http_status=r.status_code,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),receipt=str(file));s=Plain();s.feed(r.text);txt=' '.join(s.parts);(C/(name+'.txt')).write_text(txt)
  row['status']='RETRIEVED' if r.status_code==200 and len(txt)>10000 and 'Just a moment' not in txt[:100] else 'UNKNOWN_EXTERNAL_ACCESS_BLOCKED'
  row['supplement_links']=list(dict.fromkeys(link for link in s.links if ('MOESM' in link or 'sup-00' in link)))
  if name=='IDEAL_Age' and row['status']=='RETRIEVED':
   for i,link in enumerate(row['supplement_links']):
    if not link.startswith('http'):continue
    rr=requests.get(link,timeout=45);ext='.xlsx' if '.xlsx' in link else '.pdf';dest=C/(name+'_supp'+str(i+1)+ext);dest.write_bytes(rr.content)
    row.setdefault('supplement_receipts',[]).append({'url':link,'http_status':rr.status_code,'receipt':str(dest),'sha256':hashlib.sha256(rr.content).hexdigest(),'magic_hex':rr.content[:8].hex()})
 except requests.RequestException as e:row.update(status='UNKNOWN_EXTERNAL_NETWORK',error=str(e))
 return row
def main():
 with ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(fetch,SOURCES.items()))
 (O/'novelty_primary_source_receipts.json').write_text(json.dumps(rows,indent=2))
 queries=['(single-cell OR single cell OR transcriptomic) AND (aging clock OR ageing clock) AND (reliability OR test-retest OR intraclass OR reproducibility)','"sc-ImmuAging"','"RNA age metric"','"regression to the mean" AND (aging OR ageing) AND clock','clock AND responders AND rejuvenation','transcriptomic AND clock AND reliability'];pub=[]
 for i,q in enumerate(queries):
  row={'engine':'PubMed ESearch','query':q,'url':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi','queried_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
  try:
   r=requests.get(row['url'],params={'db':'pubmed','term':q,'retmode':'json','retmax':100},timeout=45);file=C/f'PubMed_{i+1}.json';file.write_bytes(r.content);j=r.json();result=j.get('esearchresult',{});row.update(http_status=r.status_code,receipt=str(file),sha256=hashlib.sha256(r.content).hexdigest(),count=result.get('count'),ids=result.get('idlist'),status='RESULTS_RETRIEVED' if result else 'UNKNOWN_EXTERNAL_API_RESPONSE')
  except (requests.RequestException,ValueError) as e:row.update(status='UNKNOWN_EXTERNAL_NETWORK',error=str(e))
  pub.append(row)
 (O/'novelty_PubMed_receipts.json').write_text(json.dumps(pub,indent=2));print(json.dumps({'primary_sources':rows,'pubmed':pub},indent=2))
if __name__=='__main__':main()
