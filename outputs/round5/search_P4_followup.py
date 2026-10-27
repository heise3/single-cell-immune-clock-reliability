#!/usr/bin/env python3
"""Sequential refinement after inspecting API errors and query semantics."""
from search_P4 import *
from urllib.parse import urlencode
import pandas as pd
def main():
 jobs={
 'HCA_projects_size50':'https://service.azul.data.humancellatlas.org/index/projects?size=50',
 'ArrayExpress_AND_longitudinal':'https://www.ebi.ac.uk/biostudies/api/v1/search?'+urlencode({'query':'PBMC AND longitudinal','collection':'arrayexpress','pageSize':100}),
 'ArrayExpress_AND_replicate':'https://www.ebi.ac.uk/biostudies/api/v1/search?'+urlencode({'query':'PBMC AND replicate','collection':'arrayexpress','pageSize':100}),
 }
 ids=json.loads((W/'GEO_longitudinal.json').read_text())['esearchresult']['idlist']
 jobs['GEO_89_summaries']='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urlencode({'db':'gds','id':','.join(ids),'retmode':'json'})
 for gse in ['GSE190992','GSE247277','GSE304704','GSE158769']:
  jobs[gse+'_SOFT']=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/soft/{gse}_family.soft.gz'
 for pmc in ['PMC10041512','PMC8983398','PMC11754665','PMC7816397','PMC8887138']:
  jobs[pmc+'_fullTextXML']=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'
 receipt=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for f in as_completed([pool.submit(get,x) for x in jobs.items()]):
   z=f.result();receipt.append(z);print(z['source'],z.get('http_status'),z.get('bytes'),z.get('error',''),flush=True)
 (O/'P4_API_followup_receipts.json').write_text(json.dumps(receipt,indent=2))
 data=json.loads((W/'cellxgene_collections.txt').read_text());rows=[]
 for d in data:
  if re.search('longitud|variation study|comparative analysis of commercial|repeat|vaccine|vaccination|PALMO|healthy.*blood|PBMC',d.get('name','')+' '+d.get('description',''),re.I):
   rows.append(d)
 (O/'P4_CELLxGENE_candidate_metadata.json').write_text(json.dumps(rows,indent=2));print('CELLXGENE_TOTAL',len(data),'candidate',len(rows),flush=True)
 for d in rows:print(d.get('collection_id'),d.get('name'),str(d.get('description'))[:180],flush=True)
if __name__=='__main__':main()
