#!/usr/bin/env python3
"""Verify fixed inputs, package targeted results, and preserve historical decisions."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import sys,json,hashlib,platform,datetime,gzip,shutil,subprocess,zipfile,importlib.metadata as meta
import numpy as np,pandas as pd,h5py,anndata as ad
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';W=P/'work/round5';R=P/'outputs/round4';A=P/'outputs/round2_A1'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 now=datetime.datetime.now(datetime.timezone.utc).isoformat();qa={};old=pd.read_csv(R/'SHA256SUMS.csv');bad=[]
 for r in old.itertuples():
  f=R/r.file
  if not f.exists() or f.stat().st_size!=r.bytes or sha(f)!=r.sha256:bad.append(r.file)
 assert not bad,bad;qa['round4_manifest_files_verified_unchanged']=len(old)
 # Actual count files used this round: two K0 halves and small subsets of AIDA.
 entries=[]
 for n in [R/'R2_K0_S_evaluation_raw_manifest.json',P/'outputs/round3/R2_E_evaluation_raw_manifest.json',P/'outputs/round3/R2_AIDA_evaluation_raw_manifest.json']:
  z=json.loads(n.read_text());assert z['complete'];entries.extend(z['chunks'])
 def verify(r):
  f=P/r['path'];d=sha(f);assert d==r['sha256'],str(f);return {'file':r['path'],'bytes':f.stat().st_size,'sha256':d,'role':'pinned_count_cache_verified_unchanged'}
 with ThreadPoolExecutor(max_workers=2) as ex:inputs=list(ex.map(verify,entries))
 qa['count_cache_files_verified_unchanged']=len(inputs);print('COUNT_HASHES_VERIFIED',len(inputs),flush=True)
 fixed=[P/'outputs/round4_checkpoint.zip',A/'main_predictions_with_LODO.csv',A/'main_predictions.csv',A/'variance_decomposition.csv',A/'compute_reliability.py',A/'run_G1_prime.py',A/'model_gene_identity_recheck.csv',P/'outputs/scImmuAging_coefficients.csv',P/'outputs/round3/AIDA_source_var.csv',P/'outputs/round3/Census_global_var.csv',P/'work/sources/GSE271896_family.soft.gz',P/'work/round2_A1/cache/hgnc_current_20261026.tsv',R/'F2_a_e_full68_predictions.csv',R/'F2_paired_ICC_bootstrap.npz',R/'Figure4_F2_paired_delta_ICC.csv',R/'F2_decision.json',W/'P4_selected_cells.parquet',W/'P4_additional_measured_counts.h5ad',W/'AIDA_cell_library_link.parquet',W/'AIDA_QC_cell_library_link.parquet']+list((R/'models').glob('*_a.npz'))+list((R/'models').glob('*_e.npz'))+list((A/'provenance_inputs').glob('K0_primary*.csv'))
 for f in fixed:inputs.append({'file':str(f.relative_to(P)),'bytes':f.stat().st_size,'sha256':sha(f),'role':'fixed_reference_or_round5_pinned_input'})
 pd.DataFrame(inputs).to_csv(O/'input_files_SHA256.csv',index=False)
 assert sha(P/'outputs/round4_checkpoint.zip')=='258655b4f4a7f81522c3567a25fbaa66f5d85b2020bc3a629b46efd322d5ff1c'
 assert len(pd.read_csv(O/'fixed_three_model_predictions.csv'))==4800;assert len(pd.read_csv(O/'P1_e_vs_published_paired_delta.csv'))==6;h=pd.read_csv(O/'P1_F2_Holm6.csv');assert h.passes_Holm6.sum()==4 and h.passes_centered_Holm6_sensitivity.sum()==2
 # Prior bootstraps were recovered with refitLODO, not substituted CIs.
 new=np.load(W/'P1_paired_bootstrap_ICCs.npz');prev=np.load(R/'F2_paired_ICC_bootstrap.npz');qa['a_e_bootstrap_max_absolute_disagreement']=max(float(abs(new[k]-prev[k]).max()) for k in prev.files);assert qa['a_e_bootstrap_max_absolute_disagreement']<1e-11
 q=pd.read_csv(O/'P2_covariate_between_donor_R2.csv');assert len(q)==15 and q.n_donors.eq(68).all();assert len(pd.read_csv(O/'P2_CMV_sex_adjusted_ICC.csv'))==30
 v=pd.read_csv(O/'P3_three_model_variance_components.csv');assert len(v)==30 and v.fit_converged.all();ff=[c+'_fraction' for c in ['donor','batch','visit_nested_plus_unmodelled_residual','cell_sampling','pseudocell_algorithm']];assert np.allclose(v[ff].sum(axis=1),1);assert v[v.model.eq('e')].pseudocell_algorithm_variance.eq(0).all();oldv=pd.read_csv(A/'variance_decomposition.csv').query('clock=="sc-ImmuAging"').sort_values(['celltype','metric']);newv=v[v.model.eq('published')].sort_values(['celltype','metric']);assert np.allclose(oldv[ff].to_numpy(),newv[ff].to_numpy(),atol=1e-12)
 e=pd.read_csv(O/'P4_external_technical_ICC.csv');assert len(e)==20 and e.n_donors.eq(20).all() and e.bootstrap_valid.eq(2000).all();assert len(pd.read_csv(O/'P4_external_technical_predictions.csv'))==400
 # Cell identity and full-library denominators; projection completeness check.
 obs=pd.read_parquet(W/'P4_selected_cells.parquet').sort_values('soma_joinid');extra=ad.read_h5ad(W/'P4_additional_measured_counts.h5ad');assert len(obs)==36805 and np.array_equal(extra.obs.soma_joinid,obs.soma_joinid);sums={}
 for r in json.loads((P/'outputs/round3/R2_AIDA_evaluation_raw_manifest.json').read_text())['chunks']:
  a=ad.read_h5ad(P/r['path'],backed='r');ix=np.flatnonzero(a.obs.soma_joinid.isin(obs.soma_joinid))
  if len(ix):
   x=a[ix].X;values=np.asarray(x.sum(axis=1)).ravel();sums.update({int(i):(float(t),float(l)) for i,t,l in zip(a.obs.soma_joinid.iloc[ix],values,a.obs.raw_sum.iloc[ix])})
  a.file.close()
 totals=np.array([sums[int(i)][0] for i in obs.soma_joinid])+np.asarray(extra.X.sum(axis=1)).ravel();lib=np.array([sums[int(i)][1] for i in obs.soma_joinid]);assert (totals<=lib).all();qa.update(P4_full_library_denominator_and_extra_gene_projection_QA=True,P1_rows=6,P2_R2_rows=15,P2_adjusted_ICC_rows=30,P3_variance_rows=30,P4_ICC_rows=20,all_bootstrap_replicates=2000,no_new_training=True,preregistration_P2_label='待声明',CD4_prior_gene_ANOVA_intermediates_disclosed=True,recorded_UTC=now)
 # Include random draws / technical cells as small metadata, never huge count matrices.
 shutil.copy2(W/'P1_paired_bootstrap_ICCs.npz',O/'P1_paired_bootstrap_ICCs.npz')
 with open(W/'P3_conditional_bootstrap_draws.csv','rb') as f,gzip.GzipFile(filename=str(O/'P3_conditional_bootstrap_draws.csv.gz'),mode='wb',mtime=0) as g:shutil.copyfileobj(f,g)
 shutil.copy2(W/'P4_selected_cells.parquet',O/'P4_selected_cell_link.parquet')
 snapshot=O/'source_snapshots';snapshot.mkdir(exist_ok=True)
 for f in sorted((W/'P4_sources').glob('*')):
  if f.is_file():
   if f.suffix=='.gz':shutil.copy2(f,snapshot/f.name)
   else:
    with f.open('rb') as src,gzip.GzipFile(filename=str(snapshot/(f.name+'.gz')),mode='wb',mtime=0) as dst:shutil.copyfileobj(src,dst)
 logs=O/'execution_logs';logs.mkdir(exist_ok=True)
 for f in W.glob('*.log'):shutil.copy2(f,logs/f.name)
 versions={'seed':20261025,'Python':sys.version,'platform':platform.platform(),'packages':{n:meta.version(n) for n in ['numpy','pandas','scipy','scikit-learn','anndata','h5py','cellxgene-census','tiledbsoma','matplotlib','statsmodels']},'R':subprocess.run(['Rscript','--version'],capture_output=True,text=True).stdout.strip(),'Census_release':'2025-11-08','AIDA_dataset_id':'c838aec3-03ef-4398-b882-0e3912abfff0','AIDA_dataset_version_id':'9deda9ad-6a71-401e-b909-5263919d85f9','published_model_code_commit':'112f8fab10c72aec6e7a7d83762a40117a1e6e24','P1_P2_shared_bootstrap':'np.random.default_rng(np.random.SeedSequence([20261025,4202]));2000 full68 multinomial weights;refitLODO','P3_MC':'A1 rng_for stableSHA256: R5_cell_bootstrap_ae|gsm|celltype and R5_a_algorithm|gsm|celltype;100 each','P4_selection':'SeedSequence([20261025,5401]);20 metadata/QC donors','P4_pseudocells':'A1 rng_for(R5_AIDA_external,gsm,celltype);100x15','P4_bootstrap':'A1 rng_for(donor_bootstrap,model,celltype);2000donor draws;refitLODO','no_new_clock_training':True,'recorded_UTC':now};(O/'versions_and_seed.json').write_text(json.dumps(versions,indent=2))
 (O/'README_reproduce.md').write_text('''# Round5 重现说明

固定种子20261025。从本项目根目录运行，Python环境为work/analysis_env/bin/python；版本见versions_and_seed.json。依赖round4_checkpoint及其绑定的A1/Round3输入与公开count缓存；input_files_SHA256.csv逐项列出输入和实际SHA，143个count缓存本轮重新校验。新checkpoint不复制约20GB矩阵，也不包含新训练。

1. 先恢复并校验round4_checkpoint；任何新结果前保留round5_execution_contract.json和gate_ledger_before_round5.json。
2. extract_GEO_covariates.py 从归档SOFT重现CMV/sex提取；analyze_P1_P2.py计算LODO、配对bootstrap、两套Holm和协变量检查。P1_paired_bootstrap_ICCs.npz保留实际2000抽分布。
3. analyze_P3.py读固定K0计数：100细胞重采样、a的100算法实现。P3_conditional_bootstrap_draws.csv.gz为实际条件抽样；无新训练。若已有work/round5/P3_noise_complete.json则复用完整条件抽样。
4. search_P4.py、search_P4_followup.py、search_P4_precise.py为有依赖的先后检索；保存所有成功及失败receipt。source_snapshots/含本轮实际检索快照；实时重跑可能随数据库更新而不同。
5. check_AIDA_repeats.py读取原始h5ad的donor/sample/library元数据；prepare_P4_technical.py逐细胞绑定并在任何split-library预测前冻结20人。保留P4_external_selection_contract.json，重现时不得根据结果重抽供者。
6. analyze_P4_technical.py使用固定已发表/e模型，不训练。旧4044基因缓存缺少但来源实测的1145基因另取整数count，保留完整library分母。网络查询曾因多列区间过慢改为稀疏整行投影。写obs时曾发生AnnData索引名冲突，recover_P4_count_write.py从已经完整写出的CSR X恢复标签，逐buffer SHA一致且保留失败文件；P4_count_write_recovery_QA.json和原日志给出证据。修复后的正常重跑已直接采用无冲突索引名，不需重现此错误。
7. summarize_P4_search.py汇总候选；plot_round5.py仅从CSV作图；build_report.py生成报告和方法结果；package_round5.py校验、更新台账和归档。

P1主要p值是近似percentile-bootstrap尾概率，不是严格随机化p。中心化零假设敏感性单独Holm，不从二者择优。P1差值CI为逐项CI。P3访视与未建模技术残差不可再分。P4为同次采血、同池、通常同lane的技术重复，不能代替90天/一年生物重测；e五类ICC点估计均低于已发表模型，不声称普遍优势得到外部复现。

预注册P2验证性/探索性标签仍为待声明，详见CD4_gene_variance_prior_exposure_declaration.md。已有CD4基因ICC及ANOVA中间量的计算与保存必须披露；本轮步骤计算完成状态不等于验证性标签。
''')
 figures=[]
 for f in sorted(O.glob('Figure*plot_data.csv')):figures.append({'file':f.name,'rows':len(pd.read_csv(f)),'sha256':sha(f)})
 pd.DataFrame(figures).to_csv(O/'plot_data_manifest.csv',index=False);(O/'final_QA.json').write_text(json.dumps(qa,indent=2))
 ledger=json.loads((P/'outputs/gate_ledger.json').read_text());ledger['round5'].update(status='COMPLETE',completed_UTC=now,seed=20261025,decision=json.loads((O/'round5_decision.json').read_text()),final_QA='outputs/round5/final_QA.json');(P/'outputs/gate_ledger.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2));shutil.copy2(P/'outputs/gate_ledger.json',O/'gate_ledger_after_round5.json')
 manifest=[]
 for f in sorted(O.rglob('*')):
  if f.is_file() and '__pycache__' not in f.parts and f.name!='SHA256SUMS.csv':manifest.append({'file':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':sha(f)})
 pd.DataFrame(manifest).to_csv(O/'SHA256SUMS.csv',index=False);archive=P/'outputs/round5_checkpoint.zip'
 with zipfile.ZipFile(archive.with_suffix('.part.zip'),'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for row in manifest:z.write(O/row['file'],arcname='round5/'+row['file'])
  z.write(O/'SHA256SUMS.csv',arcname='round5/SHA256SUMS.csv')
 archive.with_suffix('.part.zip').replace(archive)
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None
 receipt={'archive':str(archive.relative_to(P)),'bytes':archive.stat().st_size,'sha256':sha(archive),'manifest_files':len(manifest),'CRC_verified':True,'seed':20261025,'status':'COMPLETE','preregistration_P2_label':'待声明','historical_archives_unchanged':True,'recorded_UTC':now};(P/'outputs/round5_checkpoint_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2));print(json.dumps(receipt,ensure_ascii=False,indent=2),flush=True)
if __name__=='__main__':main()
