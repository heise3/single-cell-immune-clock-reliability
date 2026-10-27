#!/usr/bin/env python3
"""Public metadata searches only; no login, terms acceptance or access request."""
from pathlib import Path
import asyncio,aiohttp,json,datetime,hashlib
from urllib.parse import urlencode
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1/cache'
async def main():
 queries=[('Zenodo','GET','https://zenodo.org/api/records?'+urlencode({'q':'"siAgeModel" OR "syn61609846"','size':25}),None),('figshare','POST','https://api.figshare.com/v2/articles/search',{'search_for':'siAgeModel','limit':100}),('figshare_siAge','POST','https://api.figshare.com/v2/articles/search',{'search_for':'"siAge" "single-cell"','limit':100})];rows=[]
 async with aiohttp.ClientSession(trust_env=True,timeout=aiohttp.ClientTimeout(total=45)) as session:
  async def query(r):
   name,method,url,body=r
   try:
    async with session.request(method,url,json=body) as resp:
     text=await resp.text();(W/('siAge_'+name+'_response.txt')).write_text(text);rec={'channel':name,'url':url,'method':method,'request_body':body,'http_status':resp.status,'response_sha256':hashlib.sha256(text.encode()).hexdigest(),'reference_acquired':False}
     if resp.status==200:
      data=json.loads(text)
      hits=data.get('hits',{}).get('hits',[]) if isinstance(data,dict) else data;rec['n_hits']=len(hits);rec['hits']=[{'id':z.get('id'),'title':z.get('title') or z.get('metadata',{}).get('title'),'url':z.get('url') or z.get('links',{}).get('html')} for z in hits]
      rec['interpretation']='Review returned records; a keyword hit is not proof of native annotation-reference availability.'
     else:rec['n_hits']=None;rec['interpretation']='External access/query failure; not evidence of no deposit.'
     return rec
   except Exception as e:return {'channel':name,'url':url,'error':repr(e),'n_hits':None,'reference_acquired':False,'interpretation':'External access/query failure; not evidence of no deposit.'}
  rows=await asyncio.gather(*[query(q) for q in queries])
 report={'seed':20261025,'recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'routes':rows};(O/'siAge_alternate_API_queries.json').write_text(json.dumps(report,indent=2,ensure_ascii=False));print(json.dumps(report,indent=2,ensure_ascii=False))
asyncio.run(main())
