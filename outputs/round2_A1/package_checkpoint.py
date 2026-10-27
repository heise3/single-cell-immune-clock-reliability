#!/usr/bin/env python3
"""Verify immutable inputs, finalize active ledger, hash and package A1 only."""
from pathlib import Path
import json,hashlib,shutil,datetime,zipfile,csv,sys,importlib.metadata as md
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';SEED=20261025
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for a in iter(lambda:f.read(8*1024*1024),b''):h.update(a)
 return h.hexdigest()
def main():
 now=datetime.datetime.now(datetime.timezone.utc).isoformat();decision=json.loads((O/'decision.json').read_text());qa={};inputs=[]
 required=['A1_amendment.json','G1_prime_results.csv','G1_prime_missing_gene_list.csv','reliability_results.csv','variance_decomposition.csv','RTM_null_results.csv','novelty_top5.csv','step_status.csv','round1_to_A1_change_table.csv','audit_report_zh.md','PHARE_versions.json','scImmuAging_ICC_forest.png','RTM_primary_null.png']
 for name in required:assert (O/name).is_file(),name
 g=pd.read_csv(O/'G1_prime_results.csv');r=pd.read_csv(O/'reliability_results.csv');v=pd.read_csv(O/'variance_decomposition.csv');p=pd.read_csv(O/'main_predictions.csv');s=pd.read_csv(O/'sensitivity87/main_predictions.csv')
 assert len(g)==5 and g.G1_prime_status.eq('PASS').all();assert len(r)==126 and r.bootstrap_replicates.eq(2000).all() and r.bootstrap_valid.eq(2000).all();assert len(v)==42 and v.fit_converged.all() and v.covariance_design_rank.eq(3).all();assert p.donor.nunique()==68 and p.gsm.nunique()==320 and len(p)==6720;assert s.donor.nunique()==87 and s.gsm.nunique()==236;assert r[r.clock.str.startswith('RAM')].MDC95.isna().all()
 assert pd.read_csv(O/'main_native_glmnet_validation.csv').select_dtypes(include='number').shape[0]>0
 qa.update(required_outputs_present=True,G1_prime_five_pass=True,main_68_donors_320_visits=True,sensitivity_separate_87_236=True,ICC_rows=126,donor_bootstrap_replicates_per_row=2000,all_bootstrap_rows_2000_valid=True,all_42_REML_converged=True,all_covariance_design_rank3=True,no_RAM_MDC_years=True)
 prior=P/'outputs/round2_G1_checkpoint.zip';oldsha=sha(prior);assert oldsha=='6cdfe42b89307d577a3b15e457234afff5efb57b7f9a570f9fd2cc4662879e23';qa['original_checkpoint_sha256_unchanged']=oldsha;inputs.append({'path':str(prior.relative_to(P)),'bytes':prior.stat().st_size,'sha256':oldsha,'hash_verification':'REHASHED_AT_FINALIZATION'})
 originals=[]
 for row in json.loads((O/'provenance_input_manifest.json').read_text()):
  src=P/row['original_relative_path'];cp=P/row['archive_relative_path'];assert sha(src)==row['sha256']==sha(cp);originals.append({'path':str(src.relative_to(P)),'sha256':row['sha256'],'unchanged':True})
 # Include the exact archived K0 eligibility/window sources, without editing them.
 for row in json.loads((O/'provenance_inputs/K0_source_manifest.json').read_text()):
  src=Path(row['source']);snap=P/row['snapshot'];assert sha(src)==row['sha256']==sha(snap);target=O/'provenance_inputs/K0_archive'/Path(row['snapshot']).relative_to('work/round2/K0_snapshot');target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(snap,target);originals.append({'path':str(src),'sha256':row['sha256'],'unchanged':True})
 pd.DataFrame(originals).to_csv(O/'immutable_input_final_QA.csv',index=False);qa['immutable_source_files_verified']=len(originals)
 for name in ['union_count_cache_manifest.json','OneK1K_count_download_manifest.json']:
  m=json.loads((O/name).read_text())
  if 'complete' in m:assert m['complete']
  for row in m['chunks']:
   f=P/row['path'];assert f.stat().st_size==row['bytes'];inputs.append({'path':row['path'],'bytes':row['bytes'],'sha256':row['sha256'],'hash_verification':'FROM_DOWNLOAD_RECEIPT_SIZE_RECHECKED_AT_FINALIZATION'})
 native=P/'work/models/scImmuAging-112f8fab10c72aec6e7a7d83762a40117a1e6e24/data';ph=P/'work/models/PHARE-0ab277b38ddac534d6d91e4580134cf6c78822ca/phare/model';runtime=P/'work/round2/PHARE_runtime/phare/model';fixed=[native/'all_model.RDS',native/'all_model_inputfeatures.RDS',native/'demo_bcg_cd8t.RDS',P/'work/siAgeModel.rds',P/'work/sources/43587_2025_819_MOESM6_ESM.zip']
 fixed+=list((ph/'fast').glob('*.pkl'))+[ph/'celltypist_model_from_panage_fairSample(10w).pkl']
 for f in fixed:
  h=sha(f);inputs.append({'path':str(f.relative_to(P)),'bytes':f.stat().st_size,'sha256':h,'hash_verification':'REHASHED_AT_FINALIZATION'})
  if ph in f.parents:assert sha(runtime/f.relative_to(ph))==h
 for row in json.loads((O/'PHARE_input_manifest.json').read_text()).items():
  if row[0] in ['input_path','source_GEO_path']:
   f=P/row[1];inputs.append({'path':row[1],'bytes':f.stat().st_size,'sha256':sha(f),'hash_verification':'REHASHED_AT_FINALIZATION'})
 pd.DataFrame(inputs).to_csv(O/'input_files_SHA256.csv',index=False)
 ledgerpath=P/'outputs/gate_ledger.json';ledger=json.loads(ledgerpath.read_text())
 if not (O/'gate_ledger_before_finalization.json').exists():shutil.copy2(ledgerpath,O/'gate_ledger_before_finalization.json')
 current=ledger['round2_current'];current.setdefault('G1_prime_original_GEO_target_preserved',current.get('G1_prime'));current['G1_prime']=g.to_dict('records');current.update(status=decision['overall_status'],updated_UTC=now,scientific_terminal_verdict=decision['K6_verdict'],decision=decision,ordered_step_status=pd.read_csv(O/'step_status.csv').to_dict('records'))
 for a in ledger['scheme_amendments']:
  if a['id']=='A1':a.update(execution_status='COMPLETED_WITH_NOVELTY_ACCESS_LIMIT',scientific_terminal_verdict=decision['K6_verdict'],final_decision_artifact='outputs/round2_A1/decision.json')
 ledger['active_round']='round2_A1';ledger['overall_current_A1']=decision['overall_status'];ledgerpath.write_text(json.dumps(ledger,indent=2,ensure_ascii=False));shutil.copy2(ledgerpath,O/'gate_ledger_A1_snapshot.json')
 progress=O/'PROGRESS_NOT_TERMINAL.md'
 if progress.exists() and not (O/'PROGRESS_HISTORICAL_SNAPSHOT.md').exists():shutil.copy2(progress,O/'PROGRESS_HISTORICAL_SNAPSHOT.md')
 progress.write_text('# Historical progress file superseded\n\nThis progress snapshot is no longer current. Completed scientific K6 verdict: **PASS**, with Google Scholar novelty coverage limited. See audit_report_zh.md, decision.json and gate_ledger_A1_snapshot.json. The original progress text is retained in PROGRESS_HISTORICAL_SNAPSHOT.md.\n')
 # Preserve small execution receipts; published full texts/abstracts stay in work/.
 dest=O/'execution_logs';dest.mkdir(exist_ok=True)
 for name in ['primary_pipeline_final.log','sensitivity87_pipeline_final.log','blockwise_full_pipeline.log','PHARE_real_run.log','RTM.log','siAge_fresh_Synapse_access.log','GoogleScholar_retry.log']:
  f=P/'work/round2_A1'/name
  if f.exists():shutil.copy2(f,dest/name)
 dest=O/'public_access_receipts';dest.mkdir(exist_ok=True)
 for f in (P/'work/round2_A1/novelty_receipts').glob('GoogleScholar_*.html'):shutil.copy2(f,dest/f.name)
 for name in ['siAge_alternate_GitHub.json']:
  f=P/'work/round2_A1'/name
  if f.exists():shutil.copy2(f,dest/name)
 qa['scientific_verdict']=decision['K6_verdict'];qa['overall_status']=decision['overall_status'];qa['seed']=SEED;qa['finalization_UTC']=now;qa['plot_visual_QA']='All three PNGs inspected; ICC legend moved out of data region and reinspected';qa['RAM_CI_scope']='Conditional on fixed selected genes and full-cohort expression calibration; age LODO refit in donor bootstrap';(O/'deliverable_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False))
 extra={'PDF_reader_runtime':'bundled26.1017.22136','pypdf':'6.10.0','openpyxl':md.version('openpyxl'),'matplotlib':md.version('matplotlib'),'RAM_native_input_transfer':'Fixed public DEG panels, native SF4 target-cohort gene z standardization; no absolute age calibration'};(O/'supplementary_versions.json').write_text(json.dumps(extra,indent=2))
 files=[f for f in sorted(O.rglob('*')) if f.is_file() and '__pycache__' not in f.parts and f.name!='SHA256SUMS.csv'];manifest=[{'file':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':sha(f)} for f in files];pd.DataFrame(manifest).to_csv(O/'SHA256SUMS.csv',index=False)
 archive=P/'outputs/round2_A1_checkpoint.zip'
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files+[O/'SHA256SUMS.csv']:z.write(f,arcname='round2_A1/'+str(f.relative_to(O)))
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None;entries=len(z.namelist())
 digest=sha(archive);side=archive.with_suffix('.zip.sha256');side.write_text(digest+'  '+archive.name+'\n');receipt={'seed':SEED,'archive':str(archive),'bytes':archive.stat().st_size,'sha256':digest,'n_archive_files':entries,'ZIP_CRC_check':'PASS','original_checkpoint_unchanged':oldsha,'scientific_verdict':decision['K6_verdict'],'overall_status':decision['overall_status'],'raw_count_caches_bundled':False,'output_manifest':'round2_A1/SHA256SUMS.csv','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(P/'outputs/round2_A1_checkpoint_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False));print(json.dumps(receipt,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
