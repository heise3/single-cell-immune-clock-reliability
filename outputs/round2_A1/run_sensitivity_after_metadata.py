#!/usr/bin/env python3
from pathlib import Path
import time,subprocess
P=Path(__file__).resolve().parents[2];w=P/'work/round2_A1';deadline=time.monotonic()+7200
while not (w/'cache/HIHA_K0_cells.parquet').exists():
 if time.monotonic()>deadline:raise TimeoutError('Primary identity incomplete; no sensitivity prediction')
 time.sleep(10)
with (w/'sensitivity87_pipeline.log').open('a') as f:r=subprocess.run([str(P/'work/analysis_env/bin/python'),str(P/'outputs/round2_A1/sensitivity87/run_primary_pipeline.py')],cwd=P,stdout=f,stderr=subprocess.STDOUT)
raise SystemExit(r.returncode)
