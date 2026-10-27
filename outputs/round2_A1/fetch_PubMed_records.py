#!/usr/bin/env python3
from pathlib import Path
import json,requests,hashlib,xml.etree.ElementTree as ET
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';C=P/'work/round2_A1/novelty_receipts'
def main():
 rows=json.loads((O/'novelty_PubMed_receipts.json').read_text());warnings=[];ids=[]
 for i,row in enumerate(rows):
  j=json.loads(Path(row['receipt']).read_text()).get('esearchresult',{});warnings.append({'query':row['query'],'query_translation':j.get('querytranslation'),'warninglist':j.get('warninglist'),'errorlist':j.get('errorlist'),'returned_count':j.get('count'),'returned_ids':len(j.get('idlist',[]))})
  if i!=2:ids.extend(j.get('idlist',[])) # phrase-translation failure for RAM is not a specific negative search
 q='10.1038/s41586-025-09686-5[DOI]';r=requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'pubmed','term':q,'retmode':'json'},timeout=45);(C/'PubMed_Gong_DOI.json').write_bytes(r.content);ids.extend(r.json().get('esearchresult',{}).get('idlist',[]));ids=sorted(set(ids));r=requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',params={'db':'pubmed','id':','.join(ids),'retmode':'xml'},timeout=60);r.raise_for_status();f=C/'PubMed_records.xml';f.write_bytes(r.content);out=[]
 for a in ET.fromstring(r.content).findall('.//PubmedArticle'):
  def val(path):return ' '.join(''.join(x.itertext()) for x in a.findall(path))
  title=val('./MedlineCitation/Article/ArticleTitle');abstract=val('./MedlineCitation/Article/Abstract/AbstractText');out.append({'PMID':val('./MedlineCitation/PMID'),'title':title,'abstract':abstract,'publication_date':val('./MedlineCitation/Article/Journal/JournalIssue/PubDate'),'DOI':val('./PubmedData/ArticleIdList/ArticleId[@IdType="doi"]'),'screening_status':'TITLE_ABSTRACT_EXTRACTED_FOR_SCREENING','clock_candidate':bool('clock' in title.lower() or 'age' in title.lower())})
 full=pd.DataFrame(out);full.to_csv(C/'PubMed_title_abstract_internal.csv',index=False);public=full.drop(columns='abstract').copy();public['abstract_characters']=full.abstract.str.len();public['abstract_sha256']=full.abstract.map(lambda s:hashlib.sha256(s.encode()).hexdigest());public.to_csv(O/'novelty_PubMed_title_abstract_screen.csv',index=False);(O/'novelty_PubMed_query_warnings.json').write_text(json.dumps(warnings,indent=2));(O/'novelty_PubMed_records_receipt.json').write_text(json.dumps({'n_unique_records':len(out),'sha256':hashlib.sha256(r.content).hexdigest(),'receipt':str(f),'RAM_specific_DOI_query':q},indent=2));print(pd.DataFrame(out).loc[lambda t:t.clock_candidate,['PMID','title','DOI']].to_string(index=False))
if __name__=='__main__':main()
