#!/usr/bin/env python3
"""Finite current-round controller; never changes a gate or selects an outcome."""
from pathlib import Path
import datetime,json,time,subprocess,os,signal
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';PY=P/'work/analysis_env/bin/python'
def main():
 c=json.loads((O/'round6_execution_contract.json').read_text());deadline=datetime.datetime.fromisoformat(c['recorded_UTC'])+datetime.timedelta(minutes=60)
 while not (O/'all_count_inputs_ready.json').exists():
  if datetime.datetime.now(datetime.timezone.utc)>deadline-datetime.timedelta(seconds=120):raise TimeoutError('Frozen60minroundbudget;inputsincomplete,not abiologicalnegative')
  time.sleep(5)
 for script in ['run_predictions.py','analyze_reliability.py','plot_and_report.py']:
  remaining=(deadline-datetime.datetime.now(datetime.timezone.utc)).total_seconds();assert remaining>0
  env=os.environ.copy();env['MPLCONFIGDIR']=str(P/'work/round2/matplotlib');p=subprocess.Popen([str(PY),str(O/script)],env=env,start_new_session=True)
  try:rc=p.wait(timeout=remaining)
  except subprocess.TimeoutExpired:
   os.killpg(p.pid,signal.SIGTERM);raise TimeoutError('Frozen60minroundbudgetexhausted')
  assert rc==0,(script,rc);print('COMPLETED',script,flush=True)
if __name__=='__main__':main()
