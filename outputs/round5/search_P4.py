#!/usr/bin/env python3
"""Bounded public metadata retrieval only; no participant raw reads or clock runs."""
from pathlib import Path
import requests,json,hashlib,datetime,time,re
from concurrent.futures import ThreadPoolExecutor,as_completed
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';W=P/'work/round5/P4_sources';W.mkdir(parents=True,exist_ok=True)
JOBS={
 'cellxgene_collections':'https://api.cellxgene.cziscience.com/curation/v1/collections',
 'GEO_longitudinal':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmode=json&retmax=100&term=%28PBMC%20OR%20peripheral%20blood%20mononuclear%29%20AND%20%28single%20cell%20OR%20scRNA%29%20AND%20%28longitudinal%20OR%20repeat%20OR%20baseline%29%20AND%20%28healthy%20OR%20normal%29%20AND%20Homo%20sapiens%5BOrganism%5D%20AND%20GSE%5BETYP%5D',
 'ArrayExpress_longitudinal':'https://www.ebi.ac.uk/biostudies/api/v1/search?query=PBMC%20single%20cell%20longitudinal&pageSize=100',
 'ArrayExpress_repeat':'https://www.ebi.ac.uk/biostudies/api/v1/search?query=PBMC%20single%20cell%20replicate&pageSize=100',
 'HCA_projects':'https://service.azul.data.humancellatlas.org/index/projects?size=100',
 'GSE275071_SOFT':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275071/soft/GSE275071_family.soft.gz',
}
def get(item):
 k,url=item;t=datetime.datetime.now(datetime.timezone.utc).isoformat()
 try:
  r=requests.get(url,timeout=90,headers={'User-Agent':'Public-clock-reliability-audit/round5 (read-only)'});suffix='.gz' if url.endswith('.gz') else '.json' if ('api/' in url or 'eutils' in url or 'index/projects' in url) else '.txt';f=W/(k+suffix);f.write_bytes(r.content);return {'source':k,'query_url':url,'retrieved_UTC':t,'http_status':r.status_code,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest(),'path':str(f.relative_to(P))}
 except Exception as e:return {'source':k,'query_url':url,'retrieved_UTC':t,'error':str(e)}
def main():
 receipt=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(get,x) for x in JOBS.items()]):
   z=f.result();receipt.append(z);print(z['source'],z.get('http_status'),z.get('bytes'),z.get('error',''),flush=True)
 (O/'P4_API_retrieval_receipts.json').write_text(json.dumps(receipt,indent=2))
 f=W/'cellxgene_collections.txt'
 if f.exists():
  try:
   data=json.loads(f.read_text());hits=[{k:d.get(k) for k in ['collection_id','name','description','links','datasets','doi']} for d in data if re.search('longitud|variation study|comparative analysis of commercial|repeat|vaccine|vaccination|PALMO',str(d),re.I)];(O/'P4_CELLxGENE_candidate_metadata.json').write_text(json.dumps(hits,indent=2));print('CELLxGENE',len(data),'collections; selected',len(hits),flush=True)
  except Exception as e:print('PARSE_ERROR',repr(e),flush=True)
if __name__=='__main__':main()
