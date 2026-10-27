#!/usr/bin/env python3
"""Independent artifact/contract checks, immutable inputs, versions and ZIP."""
from pathlib import Path
import json,hashlib,datetime,shutil,zipfile,sys,platform,subprocess,importlib.metadata as md,re
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';A=P/'outputs/round2_A1';W=P/'work/round3';SEED=20261025
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 f=json.loads((O/'R2_F_decision.json').read_text());now=datetime.datetime.now(datetime.timezone.utc).isoformat();qa={};inputs=[]
 required=['F_prespecified.json','Table1_scImmuAging_reliability_wide.csv','Table2_variance_components_42.csv','R1_SEM_sensitivity_curve.csv','R1_critical_SEM_vs_CD8_estimates.csv','R2_E_reliability_results.csv','R2_accuracy_results.csv','R2_F_gate_by_type.csv','R3_overlapping_literature_2025_2026.csv','analysis_report_zh.md','manuscript_results_methods_zh.md','REPRODUCE.md','Figure_R1_SEM_tipping_point.png','Figure_R2_E_reliability.png','Figure_R2_F_gate.png']
 assert all((O/name).is_file() for name in required)
 old_manifest=pd.read_csv(A/'SHA256SUMS.csv');unchanged=[]
 for r in old_manifest.itertuples():
  h=sha(A/r.file);assert h==r.sha256,r.file;unchanged.append({'file':'outputs/round2_A1/'+r.file,'sha256':h,'unchanged':True})
 pd.DataFrame(unchanged).to_csv(O/'immutable_A1_files_final_QA.csv',index=False);archive=P/'outputs/round2_A1_checkpoint.zip';h=sha(archive);assert h=='9342f9f9562dc9c70d71daf1755e626d6f081fb05a13f6373c5dc06bdbc1ff7e';inputs.append({'path':str(archive.relative_to(P)),'bytes':archive.stat().st_size,'sha256':h,'hash_verification':'REHASHED_FINALIZATION'})
 ledgerpath=P/'outputs/gate_ledger.json';ledger=json.loads(ledgerpath.read_text());oldledger=json.loads((A/'gate_ledger_A1_snapshot.json').read_text());assert ledger['round2_current']==oldledger['round2_current'];assert ledger['scheme_amendments']==oldledger['scheme_amendments'];qa.update(A1_archive_unchanged=True,A1_manifest_files_rehashed=len(unchanged),historical_round2_and_amendments_unchanged=True)
 split=pd.read_csv(O/'R2_K0_donor_split.csv');ss=set(split.loc[split.split.eq('S'),'donor']);ee=set(split.loc[split.split.eq('E'),'donor']);assert len(ss)==len(ee)==34 and not ss&ee and len(ss|ee)==68
 short=pd.read_csv(O/'R2_S_short_pair_ledger.csv');ep=pd.read_csv(O/'R2_E_pair_ledger.csv');assert set(short.donor)<=ss and len(short)==14;assert set(ep.donor)<=ee;assert ep.groupby('stratum').size().to_dict()=={'330-400d':31,'75-105d':13}
 trsplit=pd.read_csv(O/'R2_OneK1K_donor_split.csv');train=set(trsplit.loc[trsplit.split.eq('train80'),'donor_id']);test=set(trsplit.loc[trsplit.split.eq('test20'),'donor_id']);assert len(train)==784 and len(test)==197 and not train&test
 qa.update(S_E_34_each_disjoint=True,S_short_pair_N=14,E_90d_pair_N=13,E_year_pair_N=31,OneK_split_784_197_disjoint=True)
 table=pd.read_csv(O/'Table1_scImmuAging_reliability_wide.csv');old=pd.read_csv(A/'reliability_results.csv');old=old[old.clock.eq('sc-ImmuAging')]
 assert len(table)==15
 for z in table.itertuples(index=False):
  for metric in ['prediction','age_acceleration']:
   src=old[old.celltype.eq(z.celltype)&old.stratum.eq(z.stratum)&old.metric.eq(metric)].iloc[0]
   row=table[table.celltype.eq(z.celltype)&table.stratum.eq(z.stratum)].iloc[0]
   for stat in ['ICC2_1','ICC_CI_low','ICC_CI_high','SEM','SEM_CI_low','SEM_CI_high','MDC95','MDC95_CI_low','MDC95_CI_high']:assert np.isclose(row[metric+'__'+stat],src[stat],rtol=0,atol=1e-12)
 v=pd.read_csv(O/'Table2_variance_components_42.csv');assert len(v)==42;fractions=[c for c in v if c.endswith('_fraction') and c!='largest_component_fraction'];assert np.allclose(v[fractions].sum(axis=1),1);assert v.largest_component.value_counts().to_dict()=={'batch':31,'donor':11};qa.update(R0_no_scientific_recalculation_values_verified=True,R0_table_rows=15,R0_variance_rows=42)
 curve=pd.read_csv(O/'R1_SEM_sensitivity_curve.csv');assert np.allclose(curve.SEM,np.round(np.arange(.5,6.00001,.1),1));assert curve.n_simulations.eq(100000).all();crit=curve.loc[curve[['negative_count_Holm_p','negative_median_Holm_p']].min(axis=1).lt(.05),'SEM'].max();assert crit==1.8;qa['R1_critical_grid_SEM']=crit
 r=pd.read_csv(O/'R2_E_reliability_results.csv');acc=pd.read_csv(O/'R2_accuracy_results.csv');g=pd.read_csv(O/'R2_F_gate_by_type.csv');pred=pd.read_csv(O/'R2_all_predictions.csv');assert len(r)==50 and len(acc)==75 and len(g)==15;assert r.bootstrap_replicates.eq(2000).all() and r.bootstrap_valid.eq(2000).all() and r.paired_bootstrap_valid.eq(2000).all();assert r.status.eq('ESTIMABLE').all();assert np.isfinite(pred[['prediction','age']]).all().all();assert set(pred.loc[pred.dataset.eq('Allen_E'),'donor'])==ee
 for (dataset,ct),q in pred.groupby(['dataset','celltype']):
  sets=q.groupby('variant').donor.apply(set);assert all(s==sets.iloc[0] for s in sets)
 # Re-derive the conjunction from exported underlying comparisons, not verdict text.
 independent=g.delta_ICC_vs_replicated_baseline.ge(.15)&g.MAE_increase.le(1.5)&g.variant.ne('a');assert np.array_equal(independent,g.type_meets_F);counts=g.assign(qualifies=independent).loc[g.variant.ne('a')].groupby('variant').qualifies.sum();qual=counts[counts.ge(2)].index.tolist();assert qual==f['successful_variants'];assert f['repair_verdict']==('PASS' if qual else 'FAIL');coverage=pd.read_csv(O/'R2_evaluation_feature_coverage.csv');assert coverage.nonzero_coverage.eq(1).all();qa.update(R2_models=25,R2_ICC_rows=50,R2_accuracy_rows=75,R2_all_ICC_bootstrap_valid=2000,R2_independent_F_verdict=f['repair_verdict'],R2_qualifying_variants=qual,R2_all_nonzero_gene_source_coverage=1.0)
 d=O/'models';d.mkdir(exist_ok=True)
 for ct in ['B','CD4T','CD8T','MONO','NK']:
  for var in ['a','b','c','d','e']:
   src=W/'models'/f'{ct}_{var}.npz';m=np.load(src);assert len(m['feature_ids'])==len(m['coefficients']) and np.isfinite(m['coefficients']).all() and np.isfinite(m['intercept']);dst=d/src.name;shutil.copy2(src,dst);assert sha(src)==sha(dst)
 shutil.copy2(W/'R2_ICC_bootstrap_distributions.npz',O/'R2_ICC_bootstrap_distributions.npz')
 # All raw hashes are already recorded at successful download; size is rechecked.
 for name in ['R2_OneK1K_raw_manifest.json','R2_E_evaluation_raw_manifest.json','R2_AIDA_evaluation_raw_manifest.json']:
  m=json.loads((O/name).read_text());assert m['complete']
  for z in m['chunks']:
   path=P/z['path'];assert path.stat().st_size==z['bytes'];assert z['full_row_sum_exact'];inputs.append({'path':z['path'],'bytes':z['bytes'],'sha256':z['sha256'],'hash_verification':'DOWNLOAD_HASH_SIZE_RECHECKED_FINALIZATION'})
 # Fresh hashes bind analysis aggregates, group identities, membership and solver QA inputs.
 fixed=list(W.glob('*aggregates.h5'))+list(W.glob('*groups.csv'))+list(W.glob('*cell_weights.npz'))+list((W/'training_design').glob('*.npy'))+list((W/'training_design').glob('*membership.npz'))+list(W.glob('native_Gram_*.tsv'))
 for path in fixed:inputs.append({'path':str(path.relative_to(P)),'bytes':path.stat().st_size,'sha256':sha(path),'hash_verification':'REHASHED_FINALIZATION'})
 # Prior source manifests are bound by the immutable A1 manifest rather than silently redownloaded.
 for name in ['main_Census_count_manifest.json','provenance_inputs/K0_primary_with_released_age.csv','provenance_inputs/K0_primary_pair_ledger.csv','provenance_inputs/K0_primary_batch_ledger.csv','BCG_38_paired_CD8.csv']:
  path=A/name;inputs.append({'path':str(path.relative_to(P)),'bytes':path.stat().st_size,'sha256':sha(path),'hash_verification':'REHASHED_FINALIZATION'})
 pd.DataFrame(inputs).drop_duplicates('path').to_csv(O/'input_files_SHA256.csv',index=False)
 versions=json.loads((A/'versions_and_seed.json').read_text());versions.update(python=sys.version,platform=platform.platform(),round3_seed=SEED,R1_random_stream='numpy SeedSequence([20261025,3101]) common standard normals',R2_bootstrap_stream='numpy SeedSequence([20261025,3202]) shared E donor multiplicities',other_random_streams='stable SHA256 key -> SeedSequence([20261025,first8bytes little-endian]); no Python hash() dependence',sklearn_split_random_state=SEED,bootstrap_replicates=2000,Census_release='2025-11-08',OneK_dataset_id='3faad104-2ab8-4434-816d-474d8d2641db',AIDA_dataset_id='c838aec3-03ef-4398-b882-0e3912abfff0',AIDA_version_id='9deda9ad-6a71-401e-b909-5263919d85f9',recorded_UTC=now)
 for name in ['scikit-learn','threadpoolctl','matplotlib']:versions['packages'][name]=md.version(name)
 versions['OneK_version_id']=json.loads((A/'G1_prime_execution_plan.json').read_text())['analysis_dataset_version_id']
 (O/'versions_and_seed.json').write_text(json.dumps(versions,indent=2,ensure_ascii=False));freeze=subprocess.run([sys.executable,'-m','pip','freeze'],capture_output=True,text=True,check=True);(O/'requirements_frozen.txt').write_text(freeze.stdout)
 logs=O/'execution_logs';logs.mkdir(exist_ok=True)
 for name in ['R0_R1.log','aggregate_training_v3.log','training_design.log','model_training.log','model_training_MONO.log','model_training_NK.log','model_training_CD8T.log','model_training_CD8T_LARS.log','model_training_CD4T_LARS.log','evaluation_inputs_v2.log','AIDA_source_metadata_final.log','R2_evaluation.log']:
  src=W/name
  if src.exists():shutil.copy2(src,logs/name)
 # Save historical coordinate implementation used by completed B/MONO/NK fits.
 src=W/'train_repair_models_coordinate_descent_before_acceleration.py'
 if src.exists():shutil.copy2(src,O/'train_repair_models_coordinate_descent_historical.py')
 before=O/'gate_ledger_before_finalization.json'
 if not before.exists():shutil.copy2(ledgerpath,before)
 ledger['round3_F'].update(status='ROUND3_COMPLETED',repair_verdict=f['repair_verdict'],completed_UTC=now,decision=f,ordered_step_status=pd.read_csv(O/'step_status.csv').to_dict('records'),report='outputs/round3/analysis_report_zh.md');ledger['overall_current_round3']='REPAIR_'+f['repair_verdict']+'_NOVELTY_BOUNDED_EXTERNAL_ACCURACY_LIMITED';ledgerpath.write_text(json.dumps(ledger,indent=2,ensure_ascii=False));shutil.copy2(ledgerpath,O/'gate_ledger_round3_snapshot.json')
 qa.update(seed=SEED,checked_UTC=now,scientific_status=f['repair_verdict'],plots_visual_QA='All3 PNGs inspected; interval ordering and F-target-only reference lines corrected and reinspected',large_inputs_not_bundled=True,full_published_texts_or_abstracts_not_bundled=True,access_or_numerical_failure_not_relabelled_scientific_FAIL=True)
 (O/'deliverable_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False))
 files=[p for p in sorted(O.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS.csv'];manifest=[{'file':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files];pd.DataFrame(manifest).to_csv(O/'SHA256SUMS.csv',index=False)
 # Resolve every absolute local link in the report and manuscript before delivery.
 links=[]
 for path in O.glob('*.md'):
  for target in re.findall(r'\]\((/[^)]+)\)',path.read_text()):
   assert Path(target).exists(),target;links.append(target)
 dest=P/'outputs/round3_checkpoint.zip'
 with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for path in files+[O/'SHA256SUMS.csv']:z.write(path,arcname='round3/'+str(path.relative_to(O)))
 with zipfile.ZipFile(dest) as z:assert z.testzip() is None;entries=len(z.namelist())
 digest=sha(dest);dest.with_suffix('.zip.sha256').write_text(digest+'  '+dest.name+'\n');receipt={'seed':SEED,'archive':str(dest),'sha256':digest,'bytes':dest.stat().st_size,'n_archive_files':entries,'ZIP_CRC_check':'PASS','scientific_F_verdict':f['repair_verdict'],'old_A1_archive_unchanged':h,'old_A1_files_unchanged':len(unchanged),'output_manifest':'round3/SHA256SUMS.csv','local_markdown_links_checked':len(links),'large_raw_caches_bundled':False,'recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(P/'outputs/round3_checkpoint_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False));print(json.dumps(receipt,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
