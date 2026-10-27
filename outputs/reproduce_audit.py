#!/usr/bin/env python3
"""Public-source retrieval for the ordered immune-clock feasibility audit.

No credentials, DUA acceptance, or participant-level controlled data are used.
Only gates explicitly selected on the CLI are retrieved. This is not a K4 run.
"""
import argparse, hashlib, json, re, time, urllib.parse, urllib.request, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import xml.etree.ElementTree as ET

SEED = 20261025
ROOT = Path(__file__).resolve().parent
WORK = ROOT.parent / 'work' / 'sources'
AGENT = 'ImmuneClockFeasibilityAudit/1.0 (public metadata retrieval)'
QUERIES = [
 '("single-cell aging clock"[All Fields] OR "single-cell ageing clock"[All Fields]) AND (reliability[All Fields] OR "test-retest"[All Fields] OR ICC[All Fields] OR "longitudinal stability"[All Fields])',
 '"transcriptomic clock"[All Fields] AND ("intraclass correlation"[All Fields] OR reliability[All Fields] OR "test-retest"[All Fields])',
 '"sc-ImmuAging"[All Fields] AND reproducibility[All Fields]',
 '"RNA age metric"[All Fields] AND (longitudinal[All Fields] OR "within-individual"[All Fields])',
 '("single-cell"[Title/Abstract] OR "single cell"[Title/Abstract]) AND (aging[Title/Abstract] OR ageing[Title/Abstract]) AND clock[Title/Abstract] AND (reliability[Title/Abstract] OR reproducibility[Title/Abstract] OR intraclass[Title/Abstract] OR "test-retest"[Title/Abstract])',
 '("sc-ImmuAging"[All Fields] OR "RNA Age Metric"[All Fields] OR "Tictock"[All Fields] OR "PHARE"[All Fields]) AND (aging[All Fields] OR ageing[All Fields])',
 '(single-cell[Title/Abstract] OR "single cell"[Title/Abstract]) AND (aging[Title/Abstract] OR ageing[Title/Abstract]) AND (clock[Title/Abstract] OR clocks[Title/Abstract])',
]
SCHOLAR_QUERIES=[
 '"single-cell aging clock" reliability OR "test-retest" OR ICC OR "longitudinal stability"',
 '"transcriptomic clock" "intraclass correlation"',
 '"sc-ImmuAging" reproducibility',
 '"RNA age metric" Allen longitudinal "within-individual"',
]

def fetch(url, name):
    WORK.mkdir(parents=True, exist_ok=True)
    target = WORK / name
    if target.exists():
        data = target.read_bytes()
        return {'url':url,'file':str(target.relative_to(ROOT.parent)), 'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'cached':True}
    try:
        proc=subprocess.run(['curl','-fsSL','--max-time','25','-o',str(target),url],capture_output=True,text=True,timeout=30)
        if proc.returncode: raise RuntimeError(proc.stderr.strip())
        data=target.read_bytes()
        return {'url':url,'file':str(target.relative_to(ROOT.parent)), 'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
    except Exception as exc:
        return {'url':url,'file':str(target.relative_to(ROOT.parent)),'error':str(exc)}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('gate', choices=['K2','K1'])
    args=parser.parse_args()
    records=[]
    if args.gate=='K2':
        records.append(fetch('https://www.nature.com/articles/s41586-025-09686-5','gong.html'))
        records.append(fetch('https://www.nature.com/articles/s41586-025-09686-5.pdf','gong.pdf'))
        html=(WORK/'gong.html').read_text() if (WORK/'gong.html').exists() else ''
        urls=sorted(set(re.findall(r'https://media\.springernature\.com/[^"<> ]+MOESM[^"<> ]+',html)))
        with ThreadPoolExecutor(max_workers=4) as pool:
            records.extend(pool.map(lambda url: fetch(url,urllib.parse.urlparse(url).path.rsplit('/',1)[-1]),urls))
        records.append(fetch('https://www.nature.com/articles/s43587-025-00819-z.pdf','scImmuAging.pdf'))
        records.append(fetch('https://pmc.ncbi.nlm.nih.gov/articles/PMC13403736/','pasta.html'))
        records.append(fetch('https://pubmed.ncbi.nlm.nih.gov/42525215/','biological_technical_reliability.html'))
        for i,query in enumerate(QUERIES):
            url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','term':query,'retmode':'json','retmax':200,'tool':'immune_clock_feasibility_audit'})
            record=fetch(url,f'pubmed_query_{i+1}.json')
            record['query']=query
            if 'error' not in record:
                result=json.loads((WORK/f'pubmed_query_{i+1}.json').read_text())['esearchresult']
                if int(result.get('count',0))>len(result.get('idlist',[])):
                    record=fetch(url,f'pubmed_query_{i+1}_full.json')
                    result=json.loads((WORK/f'pubmed_query_{i+1}_full.json').read_text())['esearchresult']
                    record['query']=query
                record.update(count=result.get('count'),ids=result.get('idlist'),querytranslation=result.get('querytranslation'), warnings=result.get('warninglist'))
            records.append(record); time.sleep(0.4)
        ids=set()
        for p in WORK.glob('pubmed_query_*.json'):
            ids.update(json.loads(p.read_text())['esearchresult'].get('idlist',[]))
        if ids:
            records.append(fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','id':','.join(sorted(ids)),'retmode':'xml','tool':'immune_clock_feasibility_audit'}),f'pubmed_records_{len(ids)}.xml'))
        for i,query in enumerate(SCHOLAR_QUERIES):
            record=fetch('https://scholar.google.com/scholar?'+urllib.parse.urlencode({'q':query}), 'scholar.html' if i==0 else f'scholar_{i+1}.html')
            record['query']=query;records.append(record)
    if args.gate=='K1':
        records.append(fetch('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271896&targ=self&form=text&view=full','GSE271896.txt'))
        records.append(fetch('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275067&targ=self&form=text&view=full','GSE275067.txt'))
        records.append(fetch('https://apps.allenimmunology.org/aifi/insights/dynamics-imm-health-age/','allen_insight.html'))
        records.append(fetch('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271896/soft/GSE271896_family.soft.gz','GSE271896_family.soft.gz'))
        records.append(fetch('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8387811&targ=self&form=text&view=full','GSM8387811.txt'))
    out=ROOT/f'{args.gate}_retrieval_manifest.json'
    out.write_text(json.dumps({'audit_date_local':'2026-10-25','seed':SEED,'gate':args.gate,'records':records},ensure_ascii=False,indent=2))
    for item in records:
        print(json.dumps(item,ensure_ascii=False),flush=True)

if __name__=='__main__': main()
