#!/usr/bin/env python3
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,datetime
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';py=str(P/'work/analysis_env/bin/python')
def run(script,log):
 print('BEGIN',str(script),datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
 with (W/log).open('a') as f:r=subprocess.run([py,str(script)],cwd=P,stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError(f'ENGINEERING_STOP {script}: {r.returncode}; no scientific PASS/KILL')
 print('END',str(script),flush=True)
run(O/'fetch_union_counts_blockwise.py','union_count_download.log')
with ThreadPoolExecutor(2) as ex:
 results=[ex.submit(run,O/'run_primary_pipeline.py','primary_pipeline_final.log'),ex.submit(run,O/'sensitivity87/run_primary_pipeline.py','sensitivity87_pipeline_final.log')]
 for r in results:r.result()
print('PRIMARY_AND_SENSITIVITY_RELIABILITY_COMPLETE',flush=True)
