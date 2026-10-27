#!/usr/bin/env python3
"""Anonymous, bounded access probes for the sole registered Terekhova source."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from urllib.parse import urlencode,urlsplit
import json,requests,hashlib,datetime,time
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';W=P/'work/round6';D=W/'source';D.mkdir(exist_ok=True)
PROJECT='896f377c-8e88-463e-82b0-b2a5409d6fe4'
def safeurl(url):
 q=urlsplit(url);return q.scheme+'://'+q.netloc+q.path
def fetch(name,url,download=False):
 start=time.monotonic();row={'source':name,'request_url':url,'retrieved_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'anonymous':True,'credentials_or_signature_supplied':False}
 try:
  r=requests.get(url,stream=True,timeout=(20,75),headers={'User-Agent':'clock-reliability-audit-round6/read-only'});row.update(http_status=r.status_code,redirect_chain=[{'status':x.status_code,'url_without_ephemeral_query':safeurl(x.url)} for x in r.history],final_url_without_ephemeral_query=safeurl(r.url),content_type=r.headers.get('Content-Type'),content_length=r.headers.get('Content-Length'))
  file=D/name;tmp=file.with_suffix(file.suffix+'.part');h=hashlib.sha256();n=0
  with tmp.open('wb') as f:
   for b in r.iter_content(1048576):
    f.write(b);h.update(b);n+=len(b)
    if download and n//16777216!=(n-len(b))//16777216:print('DOWNLOADED',name,n,flush=True)
    if time.monotonic()-start>600:raise TimeoutError('600sec source download cap')
  tmp.replace(file);row.update(bytes=n,sha256=h.hexdigest(),path=str(file.relative_to(P)),seconds=time.monotonic()-start)
 except Exception as e:row.update(error=repr(e),seconds=time.monotonic()-start)
 return row
def main():
 original=json.loads((P/'outputs/round5/P4_Terekhova_HCA_metadata.json').read_text());files=[]
 def walk(x):
  if isinstance(x,dict):
   if 'azul_url'in x and 'name'in x:files.append(x)
   else:
    for y in x.values():walk(y)
  elif isinstance(x,list):
   for y in x:walk(y)
 walk(original['projects'][0].get('contributedAnalyses'));ref={x['name']:x for x in files};(O/'HCA_declared_file_records.json').write_text(json.dumps(ref,indent=2))
 jobs=[('HCA_project_current.json','https://service.azul.data.humancellatlas.org/index/projects/'+PROJECT,False),('HCA_files_current.json','https://service.azul.data.humancellatlas.org/index/files?'+urlencode({'size':50,'filters':json.dumps({'project':{'is':[PROJECT]}})}),False),('b_cells.tar.gz',ref['b_cells.tar.gz']['azul_url'],True)]
 receipt=[]
 with ThreadPoolExecutor(max_workers=3) as pool:
  for fut in as_completed([pool.submit(fetch,*j) for j in jobs]):
   z=fut.result();receipt.append(z);(O/'K_a_anonymous_access_receipts.json').write_text(json.dumps(receipt,indent=2));print('ACCESS',z['source'],z.get('http_status'),z.get('bytes'),z.get('error',''),flush=True)
if __name__=='__main__':main()
