#!/usr/bin/env python3
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,subprocess,pandas as pd,numpy as np,datetime
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';W=P/'work/round6';CTS=['B','CD4T','CD8T','MONO','NK']
def main():
 assert (O/'all_count_inputs_ready.json').exists();subprocess.run([str(P/'work/analysis_env/bin/python'),str(O/'map_published_features.py')],check=True);start={'status':'RUN_STARTED_AFTER_ALL_K_AND_COUNT_QA','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'model':'published sc-ImmuAging only','no_e_c_or_training':True};(O/'prediction_run_receipt.json').write_text(json.dumps(start,indent=2))
 def one(ct):
  source=W/'extracted'/('b_cells_rna.rds' if ct=='B' else ct+'_rna.rds');target=W/f'published_{ct}_per_cell_scores.csv.gz'
  if target.exists():return
  with (W/f'prediction_{ct}.log').open('w') as f:r=subprocess.run(['Rscript',str(O/'predict_published.R'),ct,str(source),str(W/f'selected_{ct}_metadata.csv.gz'),str(O),str(W)],stdout=f,stderr=subprocess.STDOUT,timeout=900)
  assert r.returncode==0,(ct,r.returncode)
  t=pd.read_csv(O/f'native_prediction_check_{ct}.csv');cs=[c for c in t if c.startswith('cell_')];b=float((O/f'published_{ct}_intercept.txt').read_text());p=t[cs].to_numpy().T@t.coefficient.to_numpy()+b;rr=pd.read_csv(O/f'native_prediction_check_{ct}_scores.csv').prediction.to_numpy();assert np.allclose(p,rr,rtol=0,atol=1e-10),(ct,p,rr);print('PREDICT_DONE',ct,float(abs(p-rr).max()),flush=True)
 with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,CTS))
 (O/'all_predictions_complete.json').write_text(json.dumps({'status':'COMPLETE','clock':'published sc-ImmuAging only','native_linear_prediction_QA':'Python vsR doubleprecision<=1e-10','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
if __name__=='__main__':main()
