#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,datetime
import pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';W=P/'work/round2_A1';py=str(P/'work/analysis_env/bin/python')
assert pd.read_csv(O/'G1_prime_results.csv').G1_prime_status.eq('PASS').all()
steps=['fetch_HIHA_metadata.py','fetch_main_Census_counts.py','score_main_clocks.py','compute_reliability.py']
for script in steps:
 print('BEGIN_STEP',script,datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
 with (W/(script.removesuffix('.py')+'.pipeline.log')).open('a') as log:result=subprocess.run([py,str(O/script)],cwd=P,stdout=log,stderr=subprocess.STDOUT)
 if result.returncode:raise SystemExit(f'ENGINEERING_STOP {script} exit={result.returncode}; no scientific PASS/KILL inferred')
 print('END_STEP',script,flush=True)
print('PRIMARY_PIPELINE_COMPLETE',flush=True)
