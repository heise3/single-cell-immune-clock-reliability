#!/usr/bin/env python3
"""Ordinary anonymous requests only. Save access failures as failures, not no hits."""
from pathlib import Path
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
import requests,json,datetime,hashlib,re
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';C=P/'work/round2_A1/novelty_receipts';C.mkdir(exist_ok=True)
QUERIES=['"transcriptomic clock" "intraclass correlation"','"sc-ImmuAging" reproducibility','"RNA age metric" Allen longitudinal "within-individual"','"regression to the mean" aging clock','clock "responders" rejuvenation','transcriptomic clock reliability']
def query(item):
 i,q=item;url='https://scholar.google.com/scholar?'+urlencode({'hl':'en','q':q})
 row={'engine':'Google Scholar','query':q,'url':url,'queried_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  r=requests.get(url,timeout=40);name=f'GoogleScholar_{i+1}.html';(C/name).write_bytes(r.content);row.update(http_status=r.status_code,receipt=str(C/name),sha256=hashlib.sha256(r.content).hexdigest(),bytes=len(r.content))
  challenge=bool(re.search(r'unusual traffic|not a robot|captcha|automated quer|Forbidden',r.text,re.I));results=('class="gs_r gs_or gs_scl"' in r.text or 'class="gs_ri"' in r.text)
  row['status']='RESULTS_RETRIEVED' if r.status_code==200 and results and not challenge else 'UNKNOWN_EXTERNAL_ACCESS_BLOCKED';row['is_negative_search']=False;row['html_result_blocks']=r.text.count('class="gs_ri"')
 except requests.RequestException as e:row.update(status='UNKNOWN_EXTERNAL_TIMEOUT_OR_NETWORK',error=str(e),is_negative_search=False)
 return row
def main():
 with ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(query,enumerate(QUERIES)))
 (O/'GoogleScholar_retry_receipts.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False));print(json.dumps(rows,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
