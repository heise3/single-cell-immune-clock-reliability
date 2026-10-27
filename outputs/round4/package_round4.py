#!/usr/bin/env python3
"""Verify scientific export semantics, freeze receipts, preserve prior archives."""
from pathlib import Path
import json,hashlib,datetime,sys,shutil,zipfile,platform,re,importlib.metadata as md
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';W=P/'work/round4';R=P/'outputs/round3';A=P/'outputs/round2_A1';SEED=20261025
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 now=datetime.datetime.now(datetime.timezone.utc).isoformat();qa={};inputs=[];originals=[];old_archives={}
 for label,folder,digest in [('round3',R,'01622a44f2a6e6f3c38a645f6ee90a1e87989de7b0d2efa9e37234632be2240a'),('round2_A1',A,'9342f9f9562dc9c70d71daf1755e626d6f081fb05a13f6373c5dc06bdbc1ff7e')]:
  m=pd.read_csv(folder/'SHA256SUMS.csv')
  for r in m.itertuples():h=sha(folder/r.file);assert h==r.sha256,(label,r.file);originals.append({'path':str((folder/r.file).relative_to(P)),'sha256':h,'unchanged':True})
  archive=P/'outputs'/f'{label}_checkpoint.zip';assert sha(archive)==digest;old_archives[label]=digest;inputs.append({'path':str(archive.relative_to(P)),'bytes':archive.stat().st_size,'sha256':digest,'hash_verification':'REHASHED_FINALIZATION','role':'immutable recovery source'})
 pd.DataFrame(originals).to_csv(O/'immutable_prior_files_final_QA.csv',index=False);qa['prior_files_rehashed_unchanged']=len(originals);qa['old_archives_unchanged']=old_archives
 frozen=json.loads((O/'F2_prespecified.json').read_text());ledgerpath=P/'outputs/gate_ledger.json';ledger=json.loads(ledgerpath.read_text());assert all(ledger['round4_F2'][k]==v for k,v in frozen.items() if k not in ['status','repair_confirmation']);prior=json.loads((O/'gate_ledger_before_round4.json').read_text());assert all(ledger[k]==v for k,v in prior.items() if k!='active_round');assert frozen['primary_variant']=='e donor pseudobulk';assert frozen['secondary_variant']=='c genes ICC>=0.7';assert frozen['written_before_round4_calculation'];qa['F2_priority_and_historical_ledger_verified']=True
 # Assert the correction doesn't masquerade as interval-specific or pure causal decomposition.
 inventory=pd.read_csv(O/'T_31_batch_dominant_inventory.csv');assert len(inventory)==31 and inventory.clock.str.startswith('RAM').all() and inventory.interval_scope.str.contains('NOT interval-specific').all();pairs=pd.read_csv(O/'T_pair_batch_ledger.csv');counts=pairs.groupby(['stratum','batch_relation']).size().to_dict();assert counts=={('330-400d','cross_batch'):62,('7-9d','cross_batch'):15,('7-9d','same_batch'):18,('75-105d','cross_batch'):23,('75-105d','same_batch'):5};tt=pd.read_csv(O/'Figure3_batch_adjusted_reliability.csv');assert len(tt)==30 and tt.bootstrap_valid.eq(2000).all();fit=pd.read_csv(O/'T_REML_batch_components.csv');assert len(fit)==5 and fit.REML_converged.all() and fit.visits_without_external_batch_support.eq(16).all();qa.update(T31_all_RAM_pooled=True,T_same_cross_counts_verified=True,T_five_REML_converged=True,T30_reliability_rows_2000_bootstrap=True,T_pure_technical_biological_identification=False)
 # Reconstruct the primary gate using the paired bootstrap arrays themselves.
 b=np.load(W/'F2_paired_ICC_bootstrap.npz');d=pd.read_csv(O/'Figure4_F2_paired_delta_ICC.csv');decision=json.loads((O/'F2_decision.json').read_text());r=pd.read_csv(O/'Figure4_F2_full68_a_e_ICC.csv');assert len(r)==20 and r.bootstrap_valid.eq(2000).all()
 for z in d.itertuples():
  aa=b[f'{z.celltype}__a__{z.stratum}'];ee=b[f'{z.celltype}__e__{z.stratum}'];assert aa.shape==ee.shape==(2000,);lo,hi=np.percentile(ee-aa,[2.5,97.5]);assert np.allclose([lo,hi],[z.delta_CI_low,z.delta_CI_high],atol=1e-12,rtol=0);assert bool(lo>0 or hi<0)==z.CI_excludes_zero
 success=d[d.F2_target_type&d.CI_excludes_zero];assert len(success)==decision['n_target_CIs_excluding_zero']==4 and success.CI_confirms_improvement.all();assert decision['F2_result']=='确认' and decision['primary_variant']=='e';qa.update(F2_direct_paired_difference_CI_verified=True,F2_qualifying_positive_target_comparisons=4,F2_result=decision['F2_result'],F2_c_not_used=True)
 pred=pd.read_csv(O/'F2_a_e_full68_predictions.csv');assert pred.donor.nunique()==68 and pred.gsm.nunique()==320 and len(pred)==3200;identity=json.loads((O/'F2_round3_E_prediction_identity_QA.json').read_text());assert identity['matched_predictions']==1650 and identity['max_absolute_difference']<1e-10;meta=pd.read_csv(W/'K0_full_groups.csv');assert len(meta)==1600 and meta.n_cells.sum()==4609769;qa.update(full_K0_68_320=True,full_five_type_cells=4609769,old_E_fixed_predictions_recovered=True)
 # Audit every secondary fold against the donor ledger and frozen short pairs.
 splits=pd.read_csv(O/'c_50x2_donor_fold_ledger.csv');masks=pd.read_csv(O/'c_crossfit_mask_QA.csv');train=pd.read_csv(O/'c_500_model_training_QA.csv');oof=pd.read_csv(O/'c_Allen_50_repeat_OOF_predictions.csv');short=pairs[pairs.stratum.eq('7-9d')];assert len(splits)==3400 and len(masks)==len(train)==500 and len(oof)==80000;assert splits.groupby(['repeat','heldout_fold']).size().eq(34).all();assert oof.groupby(['repeat','celltype','gsm']).size().eq(1).all();assert np.isfinite(oof.prediction).all();assert train.max_dual_gap.max()<1e-5
 for rep in range(50):
  ss=splits[splits.repeat.eq(rep)]
  for held in [0,1]:
   te=set(ss.loc[ss.heldout_fold.eq(held),'donor']);tr=set(ss.loc[ss.heldout_fold.ne(held),'donor']);assert len(te)==len(tr)==34 and not te&tr;n=int(short.donor.isin(tr).sum());mm=masks[masks.repeat.eq(rep)&masks.heldout_fold.eq(held)];assert mm.n_training_short_pairs.eq(n).all();pr=oof[oof.repeat.eq(rep)&oof.heldout_fold.eq(held)];assert set(pr.donor)==te
 for ct in ['B','CD4T','CD8T','MONO','NK']:
  q=json.loads((O/f'c_{ct}_round3_recovery_QA.json').read_text());assert q['status']=='PASS' and q['max_fixed_prediction_difference_years']<1e-3
  z=np.load(W/f'c_{ct}_ICC_masks.npz');assert np.array_equal(z['mask'],z['ICC']>=.7)
 reps=pd.read_csv(O/'Figure4_c_50_repeat_ICC.csv');summary=pd.read_csv(O/'Figure4_c_repeated_crossfit_summary.csv');assert len(reps)==500 and summary.n_estimable_repeats.eq(50).all()
 for z in summary.itertuples():
  rr=reps[reps.celltype.eq(z.celltype)&reps.stratum.eq(z.stratum)].ICC2_1;assert np.allclose(np.percentile(rr,[2.5,50,97.5]),[z.ICC_repeat_p025,z.ICC_median,z.ICC_repeat_p975],atol=1e-12,rtol=0)
 qa.update(c_models=500,c_OOF_rows=80000,c_50x2_donor_separation_verified=True,c_all_repeat_ICCs_estimable=True,c_repeat_percentiles_not_bootstrap_CIs=True,c_max_dual_gap=float(train.max_dual_gap.max()))
 acc=pd.read_csv(O/'Figure5_accuracy_bias_decomposition.csv');assert len(acc)==45 and np.allclose(acc.MSE,acc.bias_squared+acc.centered_MSE,atol=1e-9);coverage=pd.read_csv(O/'AIDA_a_c_e_feature_coverage.csv');assert len(coverage)==15 and coverage.nonzero_gene_union_coverage.eq(1).all();qa.update(accuracy_rows=45,AIDA_MSE_identity_verified=True,AIDA_all_actual_nonzero_features_covered=True)
 # Bind and bundle frozen models, crossfit masks and the direct paired draws.
 modeldir=O/'models';modeldir.mkdir(exist_ok=True);cdir=modeldir/'c_crossfit';cdir.mkdir(exist_ok=True)
 for f in sorted((W/'c_models').glob('*_r*_h*')):
  if f.suffix not in ['.npz','.json']:continue
  dst=cdir/f.name;shutil.copy2(f,dst);assert sha(dst)==sha(f)
 assert len(list(cdir.glob('*.npz')))==500
 for ct in ['B','CD4T','CD8T','MONO','NK']:
  for v in ['a','e']:
   src=R/'models'/f'{ct}_{v}.npz';dst=modeldir/src.name;shutil.copy2(src,dst);assert sha(src)==sha(dst)
  shutil.copy2(W/f'c_{ct}_ICC_masks.npz',O/f'c_{ct}_ICC_masks.npz')
 shutil.copy2(W/'F2_paired_ICC_bootstrap.npz',O/'F2_paired_ICC_bootstrap.npz');shutil.copy2(R/'input_files_SHA256.csv',O/'source_input_manifest_from_round3.csv')
 # Freshly bind actual aggregate/statistic/design files read by this round.
 fixed=list(W.glob('*.h5'))+list(W.glob('*groups.csv'))+list((W/'OneK_c_stats').glob('*.npz'))+list(W.glob('*cell_weights.npz'))
 fixed +=[P/'work/round3/E_aggregates.h5',P/'work/round3/E_groups.csv',P/'work/round3/S_aggregates.h5',P/'work/round3/S_groups.csv',P/'work/round3/AIDA_aggregates.h5',P/'work/round3/AIDA_groups.csv',P/'work/sources/scImmuAging.html']
 fixed += list((P/'work/round3/training_design').glob('*pseudocells.npy'))+list((P/'work/round3/training_design').glob('*donor_means.npy'))+list((P/'work/round3/training_design').glob('*donor_metadata.csv'))
 fixed +=[A/'main_predictions_with_LODO.csv',A/'reliability_results.csv',A/'provenance_inputs/K0_primary_with_released_age.csv',A/'provenance_inputs/K0_primary_pair_ledger.csv',A/'provenance_inputs/K0_primary_batch_ledger.csv',R/'R2_training_gene_ranks.json',R/'R2_model_training_results.csv',R/'R2_reference_genes.csv',R/'R2_all_predictions.csv']
 for path in fixed:inputs.append({'path':str(path.relative_to(P)),'bytes':path.stat().st_size,'sha256':sha(path),'hash_verification':'REHASHED_FINALIZATION','role':'actual computation input or intermediate'})
 for name in ['R2_K0_S_evaluation_raw_manifest.json','E_short_additional_download_manifest.json']:
  m=json.loads((O/name).read_text());assert m['complete']
  for z in m['chunks']:
   f=P/z['path'];assert f.stat().st_size==z['bytes'];inputs.append({'path':z['path'],'bytes':z['bytes'],'sha256':z['sha256'],'hash_verification':'DOWNLOAD_OR_AGGREGATION_HASH_SIZE_RECHECKED_FINALIZATION','role':'projected count cache; full libraries preserved'})
 pd.DataFrame(inputs).drop_duplicates('path').to_csv(O/'input_files_SHA256.csv',index=False)
 versions=json.loads((R/'versions_and_seed.json').read_text());versions.update(round4_seed=SEED,python=sys.version,platform=platform.platform(),round4_primary_variant='e',round4_secondary_variant='c',round4_F2_stream='SeedSequence([20261025,4202]); shared full68 donor multiplicities',round4_c_split_stream='StratifiedKFold2 shuffle; random_state=SeedSequence([20261025,4401,repeat]).generate_state(1)[0]',round4_pseudocell_stream='Reuse exact E membership; remaining cache uses stable SHA256 keyed SeedSequence as frozen round3 native preprocessing',T_bootstrap_stream='A1 stable SHA256 rng_for(clock=sc-ImmuAging_crossdonor_batch_adjusted,celltype);2000',T_CI_scope='conditional on fitted cross-donor batch correction',c_uncertainty='50 repeated splits, not donor-bootstrap CI',recorded_UTC=now)
 for name in ['numpy','pandas','scipy','scikit-learn','h5py','anndata','cellxgene-census','tiledbsoma','matplotlib','threadpoolctl']:versions['packages'][name]=md.version(name)
 (O/'versions_and_seed.json').write_text(json.dumps(versions,indent=2,ensure_ascii=False));shutil.copy2(R/'requirements_frozen.txt',O/'requirements_frozen.txt')
 logs=O/'execution_logs';logs.mkdir(exist_ok=True)
 for f in W.glob('*.log'):
  if f.name.startswith('package'):continue
  shutil.copy2(f,logs/f.name)
 for f in logs.glob('package*.log'):f.unlink()
 (logs/'README.md').write_text('''# Execution log status

Accepted completion logs: T_analysis.log, prepare_extra_contiguous.log, prepare_c_stats_final.log, c_masks.log, c_fit_*.log, S_prefetch0/1.log, confirmation_pipeline.log, plots_delivery.log. The checkpoint receipt and deliverable_QA.json verify packaging.

Historical engineering logs are retained: prepare_S.log / prepare_extra.log show the initial sandbox proxy denial, before the authorized public-data network run. prepare_S_network.log and prepare_extra_network.log were interrupted to switch to disjoint download workers / full-row reads with identical feature projection. prepare_c_stats.log records a metadata column-name error fixed before accepted fold/statistic outputs; prepare_c_stats_final.log is the successful rerun. These are not negative scientific results and do not supersede completion markers. No gate, cohort, model priority, threshold or outcome was changed in response.
''')
 plots=[]
 for n in [2,3,4,5]:
  for f in sorted(O.glob(f'Figure{n}*.csv')):plots.append({'figure':n,'data_file':f.name,'n_rows':len(pd.read_csv(f)),'scope':'See figure_captions_zh.md for source/statistic definitions'})
 pd.DataFrame(plots).to_csv(O/'plot_data_manifest.csv',index=False);chrono=[]
 for name in ['F2_prespecified.json','round4_execution_contract.json','c_50x2_donor_fold_ledger.csv','c_crossfit_mask_QA.csv','F2_decision.json','c_500_model_training_QA.csv']:
  f=O/name;chrono.append({'file':name,'sha256':sha(f),'mtime_UTC':datetime.datetime.fromtimestamp(f.stat().st_mtime,datetime.timezone.utc).isoformat()})
 (O/'chronology_receipt.csv').write_text(pd.DataFrame(chrono).to_csv(index=False));ledger['round4_F2'].update(status='ROUND4_COMPLETED',repair_confirmation=decision['F2_result'],decision=decision,ordered_step_status=pd.read_csv(O/'step_status.csv').to_dict('records'),completed_UTC=now,report='outputs/round4/analysis_report_zh.md');ledger['overall_current_round4']='F2_CONFIRMED_INTERNAL_EXTENSION_TECHNICAL_IDENTIFICATION_LIMITED';ledgerpath.write_text(json.dumps(ledger,indent=2,ensure_ascii=False));shutil.copy2(ledgerpath,O/'gate_ledger_round4_snapshot.json')
 qa.update(seed=SEED,checked_UTC=now,all_four_figures_visual_QA='Inspected; legends moved outside data region and figures4/5 rechecked',large_count_caches_bundled=False,models_bundled=510,scientific_verdict=decision['F2_result']);(O/'deliverable_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False))
 for path in O.glob('*.md'):
  for target in re.findall(r'\]\((/[^)]+)\)',path.read_text()):assert Path(target).exists(),target
 files=[f for f in sorted(O.rglob('*')) if f.is_file() and '__pycache__' not in f.parts and f.name!='SHA256SUMS.csv'];manifest=[{'file':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':sha(f)} for f in files];pd.DataFrame(manifest).to_csv(O/'SHA256SUMS.csv',index=False);archive=P/'outputs/round4_checkpoint.zip'
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files+[O/'SHA256SUMS.csv']:z.write(f,arcname='round4/'+str(f.relative_to(O)))
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None;entries=len(z.namelist())
 digest=sha(archive);archive.with_suffix('.zip.sha256').write_text(digest+'  '+archive.name+'\n');receipt={'seed':SEED,'archive':str(archive),'bytes':archive.stat().st_size,'sha256':digest,'n_archive_files':entries,'ZIP_CRC_check':'PASS','F2_result':decision['F2_result'],'scope':'internal full68 extension; score-level batch adjustment, not causal decomposition','old_archives_unchanged':old_archives,'old_files_verified':len(originals),'models_bundled':510,'large_raw_count_caches_bundled':False,'manifest':'round4/SHA256SUMS.csv','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(P/'outputs/round4_checkpoint_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False));print(json.dumps(receipt,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
