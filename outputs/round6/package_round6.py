#!/usr/bin/env python3
"""Read-only validation, source binding, provenance and checkpoint packaging."""
from pathlib import Path
import json,hashlib,datetime,subprocess,sys,shutil,zipfile,tarfile,importlib.metadata
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';W=P/'work/round6'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(name,x):(O/name).write_text(json.dumps(x,ensure_ascii=False,indent=2))
def main():
 r=pd.read_csv(O/'replication_reliability_results.csv');p=pd.read_csv(O/'published_predictions_with_LODO.csv');d=json.loads((O/'replication_decision.json').read_text());pop=json.loads((O/'primary_population_freeze.json').read_text());assert len(r)==10 and (r.bootstrap_valid==2000).all() and (r.n_donors==101).all();assert len(p)==1250 and p.donor.nunique()==101 and p.groupby('celltype').gsm.nunique().eq(250).all();assert pop['n_cells']==1421481
 assert all(json.loads((O/f'K_{k}_status.json').read_text())['status'].startswith('PASS') for k in ['a','b','c'])
 qa=pd.read_csv(O/'all_type_count_reconstruction_QA.csv');assert len(qa)==5 and qa.max_reconstruction_integer_error.lt(1e-7).all() and qa.max_feature_count_vs_metadata_error.eq(0).all()
 c=pd.read_csv(O/'replication_comparison_with_AIDA.csv').set_index('celltype');assert d['replication_verdict']==('复现' if c.loc['B','ICC2_1']<=.6 and c.loc['MONO','ICC2_1']<=.6 and c.strictly_below_AIDA_technical.sum()>=3 else '未复现')
 # Independent scalar LODO check using saved multiplicities, without new resampling.
 sys.path.insert(0,str(P/'outputs/round2_A1'));from compute_reliability import lodo,batched_lodo
 draws=np.load(W/'donor_bootstrap_multiplicities.npy');lqa=[]
 for ct,g in p.groupby('celltype'):
  g=g.sort_values(['donor','gsm']).reset_index(drop=True);y=g.prediction.to_numpy();v=batched_lodo(g,y,draws[:3])
  for i in range(3):
   er=float(np.max(abs(v[i]-lodo(g,y,draws[i]))));assert er<1e-10;lqa.append({'celltype':ct,'draw':i,'max_abs_error':er})
 write('LODO_refit_numerical_QA.json',{'status':'PASS','checks':lqa,'no_new_resampling':True})
 # Pin reused artifacts to round5, and hash all actual input matrices.
 refs=['outputs/round5/P4_external_technical_ICC.csv','outputs/round5/P2_chronological_age_signal.csv','outputs/round2_A1/compute_reliability.py','outputs/round2_A1/run_G1_prime.py','outputs/round2_A1/model_gene_identity_recheck.csv','outputs/scImmuAging_coefficients.csv','work/round2_A1/cache/hgnc_current_20261026.tsv','outputs/round2_A1/published_model_gene_annotation.csv','outputs/round2_A1/OneK1K_var.csv']
 old=pd.read_csv(P/'outputs/round5/input_files_SHA256.csv').set_index('file');sums=pd.read_csv(P/'outputs/round5/SHA256SUMS.csv').set_index('file');checks=[]
 for f in refs:
  actual=sha(P/f);expected=None
  if f in old.index:expected=old.loc[f,'sha256']
  elif f.startswith('outputs/round5/') and f.split('/')[-1] in sums.index:expected=sums.loc[f.split('/')[-1],'sha256']
  if expected is not None:assert actual==expected,f
  checks.append({'file':f,'sha256':actual,'prior_pinned_hash_verified':expected is not None})
 write('reused_reference_integrity.json',{'status':'PASS','files':checks,'unlisted_annotation_inputs_hashed_currently':True})
 manifests=json.loads((O/'primary_matrix_download_manifest.json').read_text());sources=[];binding=[]
 for m in manifests:
  assert m['status']=='COMPLETE';sources += [(m['path'],m['expected_sha256'],'anonymous_source_archive'),(m['RNA_RDS_path'],m['RNA_RDS_sha256'],'RNA_RDS')]
 sources += [('work/round6/source/b_cells.tar.gz',None,'anonymous_source_archive'),('work/round6/extracted/b_cells_rna.rds',None,'RNA_RDS')]
 inputs={f:(expected,role) for f,expected,role in sources}
 for f in refs:inputs[f]=(None,'fixed_reference_or_existing_CSV')
 for f in sorted((W/'extracted').glob('*metadata.csv')):inputs[str(f.relative_to(P))]=(None,'per_cell_metadata')
 inputs[str((W/'primary_selected_cells.parquet').relative_to(P))]=(None,'frozen_selected_cells')
 for f in sorted((W/'source').glob('*.json')):inputs[str(f.relative_to(P))]=(None,'source_receipt')
 for f in sorted((W/'source').glob('author_*')):inputs[str(f.relative_to(P))]=(None,'author_code_not_executed_or_bundled')
 for f in sorted((W/'source').glob('*.pdf')):inputs[str(f.relative_to(P))]=(None,'annual_design_source_not_bundled')
 def hashed(item):
  f,(expect,role)=item;path=P/f;actual=sha(path)
  if expect:assert actual==expect,f
  return {'file':f,'bytes':path.stat().st_size,'sha256':actual,'role':role,'declared_hash_verified':bool(expect)}
 with ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(hashed,inputs.items()))
 # Bind metadata obtained by prefix streaming to the final SHA-verified archives.
 for ct,archive,meta in [('B',W/'source/b_cells.tar.gz',W/'extracted/b_cells_metadata.csv')]+[(m['celltype'],P/m['path'],W/f"extracted/{m['celltype']}_metadata.csv") for m in manifests]:
  with tarfile.open(archive,'r|gz') as tf:
   for member in tf:
    if member.name.endswith('_metadata.csv'):
     h=hashlib.sha256();f=tf.extractfile(member)
     for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
     digest=h.hexdigest();assert digest==sha(meta),ct;binding.append({'celltype':ct,'archive_member':member.name,'metadata_sha256':digest,'matched':True});break
   else:raise RuntimeError('metadata member absent')
 write('metadata_archive_binding_QA.json',{'status':'PASS','files':binding});pd.DataFrame(rows).to_csv(O/'input_files_SHA256.csv',index=False)
 shutil.copy2(W/'donor_bootstrap_multiplicities.npy',O/'donor_bootstrap_multiplicities.npy')
 logdir=O/'execution_logs';logdir.mkdir(exist_ok=True)
 for f in W.glob('*.log'):shutil.copy2(f,logdir/f.name)
 versions={'seed':20261025,'Python':sys.version,'platform':sys.platform,'packages':{x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','scikit-learn','anndata','h5py','matplotlib','statsmodels']},'R':subprocess.check_output(['Rscript','-e','cat(R.version.string,"\nMatrix ",as.character(packageVersion("Matrix")),"\n");sessionInfo()'],text=True),'pseudocell_seed_stream':'SHA256 rng_for(seed, R6_Terekhova_published,donor,visit,type),native_indices100x15','bootstrap_seed_stream':'numpy SeedSequence([20261025,6201]),shared multinomial(101) x 2000'};write('versions_and_seed.json',versions)
 now=datetime.datetime.now(datetime.timezone.utc);contract=json.loads((O/'round6_execution_contract.json').read_text());start=datetime.datetime.fromisoformat(contract['recorded_UTC']);elapsed=(now-start).total_seconds();decisiontime=datetime.datetime.fromisoformat(d['recorded_UTC']);assert (decisiontime-start).total_seconds()<=3600
 scope='本轮取得的逐细胞CSV不含日期或相对天数，实际采样日间隔未取得；按年度设计和Tube编码推断名义波次，未穷尽所有临床文件。'
 ledger=P/'outputs/gate_ledger.json';led=json.loads(ledger.read_text());before=json.loads((O/'gate_ledger_before_round6.json').read_text());unchanged={k:led.get(k)==v for k,v in before.items() if k not in ['active_round','round6']};assert all(unchanged.values());write('prior_round_ledger_integrity.json',{'status':'PASS','unchanged_keys':unchanged})
 led['round6']['closeout']={'status':'COMPLETE','recorded_UTC':now.isoformat(),'elapsed_seconds_before_archive':elapsed,'analysis_completed_seconds':(decisiontime-start).total_seconds(),'within_60min_analysis_bound':True,'metadata_time_scope_clarification':scope,'all_five_count_QA_pass':True,'same_batch_pairs':4,'cross_batch_pairs':97,'no_new_variants_or_training':True,'preregistration_P2_CD4_label':'待声明'};ledger.write_text(json.dumps(led,ensure_ascii=False,indent=2));shutil.copy2(ledger,O/'gate_ledger_after_round6.json')
 write('final_QA.json',{'status':'PASS','n_donors':101,'n_visits':250,'n_predictions':1250,'n_reliability_rows':10,'bootstrap_valid_each':2000,'metadata_time_scope':scope,'analysis_elapsed_seconds':(decisiontime-start).total_seconds(),'closeout_elapsed_seconds':elapsed,'all_preround_ledger_entries_unchanged':True,'numerical_and_input_integrity_checks_pass':True,'replication_verdict':d['replication_verdict']})
 (O/'README_REPRODUCE.md').write_text('''# Round6 reproducibility receipt\n\nSeed:20261025. Only published sc-ImmuAging; no e/c execution or training. Gate and rule contract was written before any download or data computation.\n\nRun from project root with work/analysis_env/bin/python. These scripts intentionally reference the inherited A1 helper scripts and frozen input annotations listed in input_files_SHA256.csv. Matrices remain in work/round6 and are NOT included in this small checkpoint. Recover them from anonymous HCA object URLs in HCA_declared_file_records.json / primary_matrix_download_manifest.json and verify SHA256 before use. Source README/author code/PDF are linked and hashed, not redistributed.\n\nOrder: check_Ka_access.py, inspect_B_input_format.R, complete_Ka.py; check_Kb.py; fetch_Kc_metadata.py and check_Kc.py; download_primary_matrices.py; run_count_QA.py; run_predictions.py; analyze_reliability.py; read_existing_supplement_and_batches.py; plot_and_report.py; package_round6.py. Stop on first K failure. Reuse the frozen donor/visit/cell ledgers; do not redefine cohort after predictions. Full source matrix QA preceded all predictions. The finite controller run_when_inputs_ready.py enforced the one-hour analysis deadline.\n\nNative QC is 200..5000 features, mito<=25%; all five types>=30 cells. Each donor contributes one earliest/latest healthy visit pair. All 250 complete visits of the paired 101 donors enter equal-donor-weight LODO. No single-visit-only donors enter calibration. Repeat LODO for each of 2000 shared donor bootstrap multiplicities, excluding all copies of held original donor. Raw and residual ICC(2,1), SEM, MDC95 and percentile95%CI are in replication_reliability_results.csv.\n\nTime scope: retrieved per-cell CSV lacks dates/day offsets. Annual study design plus source Tube coding supports nominal waves; exact intervals remain unknown. 39 nominal1-year and62 nominal2-year pairs are pooled by the frozen rule. Do not claim exact365/730 days or exhaustive clinical-file absence.97 cross-batch and4 same-batch pairs mean this does not isolate biological variation.\n\nAllen Pearson values were read only from the existing round5 CSV. Previously disclosed CD4 gene-level ICC/MSR/MSC/MSE exposure remains; preregistration P2 label stays 待声明. No new CD4 gene variance decomposition was run in round6.\n\nSHA256SUMS.csv lists checkpoint members excluding itself; input_files_SHA256.csv hashes external/source inputs. The sibling round6_checkpoint_receipt.json identifies the ZIP SHA and verifies CRC and every manifest member.\n''')
 files=[f for f in sorted(O.rglob('*')) if f.is_file() and '__pycache__' not in f.parts and f.name!='SHA256SUMS.csv'];pd.DataFrame([{'file':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':sha(f)} for f in files]).to_csv(O/'SHA256SUMS.csv',index=False)
 archive=P/'outputs/round6_checkpoint.zip'
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files+[O/'SHA256SUMS.csv']:z.write(f,'round6/'+str(f.relative_to(O)))
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  for row in pd.read_csv(O/'SHA256SUMS.csv').itertuples():assert hashlib.sha256(z.read('round6/'+row.file)).hexdigest()==row.sha256
 receipt={'path':str(archive),'sha256':sha(archive),'bytes':archive.stat().st_size,'manifest_files':len(files),'ZIP_CRC_and_each_manifest_SHA256':'PASS','seed':20261025,'recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(P/'outputs/round6_checkpoint_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2));print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
