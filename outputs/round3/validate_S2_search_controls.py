#!/usr/bin/env python3
from pathlib import Path
import requests,json,time,hashlib,datetime
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';D=P/'work/round3/novelty_receipts'
def main():
 out=[]
 for i,q in enumerate(['TranslAGE','"Single-cell immune aging clocks"','"transcriptomic clock"']):
  row={'query':q,'purpose':'positive control for bulk parser/index; not a no-overlap test','queried_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
  try:
   r=requests.get('https://api.semanticscholar.org/graph/v1/paper/search/bulk',params={'query':q,'year':'2025-2026','fields':'title,year,publicationDate,externalIds,url'},timeout=45);f=D/f'S2_positive_control_{i+1}.json';f.write_bytes(r.content);j=r.json();row.update(http_status=r.status_code,total_hits=j.get('total'),first_titles=[a['title'] for a in j.get('data',[])[:8]],sha256=hashlib.sha256(r.content).hexdigest(),receipt=str(f))
  except (requests.RequestException,ValueError) as e:row['error']=str(e)
  out.append(row);time.sleep(3)
 (O/'R3_S2_positive_controls.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
