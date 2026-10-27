#!/usr/bin/env python3
"""Blood-filtered HCA and exact-AND ArrayExpress refinements; inspect primary metadata."""
from search_P4 import *
from urllib.parse import urlencode
def main():
 jobs={
 'HCA_blood_p1':'https://service.azul.data.humancellatlas.org/index/projects?'+urlencode({'size':50,'filters':json.dumps({'organ':{'is':['blood']}})}),
 'ArrayExpress_PBMC_single_healthy':'https://www.ebi.ac.uk/biostudies/api/v1/search?'+urlencode({'query':'PBMC AND "single cell" AND healthy','collection':'arrayexpress','pageSize':100}),
 'ArrayExpress_PBMC_single_repeated':'https://www.ebi.ac.uk/biostudies/api/v1/search?'+urlencode({'query':'PBMC AND "single cell" AND (replicate OR longitudinal OR baseline)','collection':'arrayexpress','pageSize':100}),
 }
 for gse in ['GSE206265','GSE161329']:
  jobs[gse+'_SOFT']=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz'
 receipt=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(get,x) for x in jobs.items()]):
   z=f.result();receipt.append(z);print(z['source'],z.get('http_status'),z.get('bytes'),z.get('error',''),flush=True)
 z=json.loads((W/'HCA_blood_p1.json').read_text());nxt=z.get('pagination',{}).get('next')
 if nxt:receipt.append(get(('HCA_blood_p2',nxt)))
 (O/'P4_API_precise_receipts.json').write_text(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
