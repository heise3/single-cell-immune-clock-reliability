#!/usr/bin/env python3
"""Queue count-only checks as each pinned source object completes; no predictions."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import pandas as pd,json,time,subprocess,datetime
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';W=P/'work/round6';CTS=['B','CD4T','CD8T','MONO','NK']
def main():
 assert all(json.loads((O/f'K_{g}_status.json').read_text())['status']=='PASS' for g in ['a','b','c']);z=pd.read_parquet(W/'primary_selected_cells.parquet')
 for ct,g in z.groupby('celltype'):g.to_csv(W/f'selected_{ct}_metadata.csv.gz',index=False,compression='gzip')
 start=time.monotonic()
 def one(ct):
  source=W/'extracted'/('b_cells_rna.rds' if ct=='B' else ct+'_rna.rds')
  while ct!='B':
   man=O/'primary_matrix_download_manifest.json';d=json.loads(man.read_text()) if man.exists() else [];row=next((r for r in d if r['celltype']==ct),None)
   if row:
    if row['status']=='FAILED':raise RuntimeError(row)
    if row['status']=='COMPLETE':break
   if time.monotonic()-start>2000:raise TimeoutError('Source wait budget')
   time.sleep(8)
  out=O/f'count_reconstruction_QA_{ct}.csv'
  if out.exists():return ct
  with (W/f'count_QA_{ct}.log').open('w') as f:r=subprocess.run(['Rscript',str(O/'validate_count_reconstruction.R'),ct,str(source),str(W/f'selected_{ct}_metadata.csv.gz'),str(O)],stdout=f,stderr=subprocess.STDOUT,timeout=900)
  assert r.returncode==0,(ct,r.returncode);assert out.exists();print('QA_FINISHED',ct,flush=True);return ct
 with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,CTS))
 q=pd.concat([pd.read_csv(O/f'count_reconstruction_QA_{ct}.csv') for ct in CTS],ignore_index=True);assert q.native_QC_verified.all() and not q.clock_values_computed.any();q.to_csv(O/'all_type_count_reconstruction_QA.csv',index=False);(O/'all_count_inputs_ready.json').write_text(json.dumps({'status':'PASS_ALL_FIVE_BEFORE_ANY_CLOCK_PREDICTION','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_cells':int(q.analysis_cells.sum()),'verification':'All analysis-cell releasednonzeros recoverintegercounts;full library sums and nativegene counts independentlymatch releasedQCmetadata;no scale assumptions used without this check'},indent=2));print('ALL_COUNT_INPUTS_READY',flush=True)
if __name__=='__main__':main()
