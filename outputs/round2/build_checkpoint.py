#!/usr/bin/env python3
"""Build an honest, hashed checkpoint; never fill unrun scientific results."""
from pathlib import Path
import json,hashlib,platform,importlib.metadata as md,shutil,datetime
import pandas as pd

OUT=Path(__file__).resolve().parent
PROJECT=OUT.parents[1]
WORK=PROJECT/'work/round2'
SEED=20261025
checkpoint_status={'status':'PAUSED_AT_G1_BY_USER','decision':'保持原门槛，停留 G1','scientific_terminal_verdict':'NOT_EVALUABLE',
 'terminal_KILL_triggered':False,'terminal_PASS_declared':False,'downstream_analysis_authorized_under_changed_gate':False,
 'automatic_resume':False,'updated_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(OUT/'checkpoint_status.json').write_text(json.dumps(checkpoint_status,indent=2,ensure_ascii=False))
pop=json.loads((OUT/'step0_summary.json').read_text())
g=pd.read_csv(OUT/'G1_gene_coverage_GEO.csv')
mask=json.loads((OUT/'G1_official_CD8_demo_mask_summary.json').read_text())
mapping=pd.read_csv(OUT/'scImmuAging_gene_mapping.csv').fillna('')
coef=pd.read_csv(PROJECT/'outputs/scImmuAging_coefficients.csv')
nonzero=coef[coef.coefficient.ne(0)&coef.gene.ne('(Intercept)')].merge(mapping,left_on='gene',right_on='model_ensembl')
nonzero[nonzero.mapped_geo_ids.eq('')].to_csv(OUT/'G1_residual_missing_gene_weights.csv',index=False)

coverage=[]
sources=[('G1_glaucoma_candidate_var.csv','GSE268936','30c2a6fd-d547-460f-a5e7-44c62a2af7ad','0be8ebd5-bb46-4bc7-a998-f24ce59bb28a'),
 ('G1_multitissue2025_candidate_var.csv','GSE299043 ILC partition','e3da8320-0b3c-4bae-ad55-a89eb59ca870','644ba356-40ca-4c75-a491-a8c1b3e47a11')]
for file,source,dataset,version in sources:
 v=pd.read_csv(OUT/file).fillna('');ids=set(v['_index'])
 for ct,z in nonzero[nonzero.mapped_geo_ids.eq('')].groupby('celltype'):
  # Report only proven exact-ID coverage. This is a feature-universe screen,
  # never a matched-lineage experiment or a gene's expression measurement.
  coverage.append({'source':source,'dataset_id':dataset,'dataset_version_id':version,'screened_feature_count':len(v),
   'celltype':ct,'residual_genes_to_mask':len(z),'residual_ids_present':int(z.gene.isin(ids).sum()),
   'residual_ids_absent':'|'.join(z.loc[~z.gene.isin(ids),'gene']),
   'status':'ANNOTATION_SCREEN_ONLY_PARTIAL_COVERAGE_NOT_G1_PASS'})
pd.DataFrame(coverage).to_csv(OUT/'G1_independent_validation_coverage.csv',index=False)

gates=[{'step':'0','status':'PASS','evidence':'68 donors, 320 visits; primary batches=90; paired donors 33/28/62 for 7–9/75–105/330–400 days',
 'conclusion':'Frozen K0 main population restored; 87-person definition remains sensitivity only.'},
 {'step':'1:scImmuAging_gene_mapping','status':'PASS_COVERAGE_COMPONENT_ONLY','evidence':'; '.join(f'{r.celltype}: {r.missing_after_Ensembl_HGNC}/{r.nonzero_genes}' for r in g.itertuples()),
 'conclusion':'Residual feature fractions below 5%; this is not proof of masking robustness.'},
 {'step':'1:G1_CD8','status':'EXCLUDE_FROM_MAIN','evidence':f"n=5, r={mask['Pearson_r']:.9f}, bias={mask['bias_years']:.6f}, delta_MAE={mask['MAE_change']:.6f}",
 'conclusion':'User confirmed all CD8 reliability/RTM results supplementary only and ineligible to support terminal KILL.'},
 {'step':'1:G1_other_types','status':'HOLD_AT_G1_BY_USER','evidence':'Inspected public independent PBMC inputs do not jointly represent the complete residual mask; no fabricated r=1 pass.',
 'conclusion':'Complete matched-type masking experiments have not been run.'},
 {'step':'2','status':'NOT_RUN','evidence':'G1 unresolved','conclusion':'No main-cohort ICC, SEM, MDC, bootstrap CI, LODO residuals, or variance decomposition computed.'},
 {'step':'3','status':'NOT_RUN_COMPLETE_COVERAGE','evidence':'PHARE model-only smoke works; siAge anonymous reference download returns 403.',
 'conclusion':'At least two fully preprocessed clocks have not completed the main analysis.'},
 {'step':'4','status':'NOT_RUN','evidence':'No eligible CD8 main SEM; CD8 supplementary amendment recorded.',
 'conclusion':'No new RTM simulations, p values, correction, or IFN-adjusted result.'},
 {'step':'5','status':'NOT_RUN','evidence':'Sequential workflow has not reached novelty update.',
 'conclusion':'No round2 claim of no complete novelty overlap.'},
 {'step':'6','status':'NOT_EVALUABLE','evidence':'Required G1/main reliability evidence incomplete.',
 'conclusion':'Neither scientific PASS nor terminal KILL declared.'}]
pd.DataFrame(gates).to_csv(OUT/'step_status.csv',index=False)
rows=[]
for r in g.itertuples():
 rows.append({'metric':'missing_nonzero_genes','celltype':r.celltype,'round1':r.missing_round1,'round2':r.missing_after_Ensembl_HGNC,
 'change':r.missing_after_Ensembl_HGNC-r.missing_round1,'interpretation':'Annotation repair only; no new main-cohort reliability estimate'})
rows.extend([{'metric':'candidate_donors','celltype':'','round1':'89 / 87','round2':68,'change':'definitions differ','interpretation':'Frozen K0 68 main; 87 sensitivity only'},
 {'metric':'main_visits','celltype':'','round1':'40 pilot visits from 20 donors','round2':'320 frozen eligible visits','change':'population alignment','interpretation':'Eligibility ledger, not completed numerical analysis'},
 {'metric':'G1_CD8_mask_r','celltype':'CD8T','round1':'NOT_RUN','round2':mask['Pearson_r'],'change':'new check','interpretation':'Excluded from main'},
 {'metric':'ICC/SEM/MDC','celltype':'all','round1':'See preserved round1 K4 tables','round2':'NOT_RUN','change':'not estimable','interpretation':'Do not recycle pilot numbers into the 68-donor main cohort'}])
pd.DataFrame(rows).to_csv(OUT/'round1_vs_round2_changes.csv',index=False)

versions={'seed':SEED,'python':platform.python_version(),'packages':{}}
for name in ['numpy','scipy','pandas','h5py','fsspec','aiohttp','requests','scanpy','celltypist','statsmodels']:
 try:versions['packages'][name]=md.version(name)
 except md.PackageNotFoundError:versions['packages'][name]=None
versions['pinned_models']={'scImmuAging':'112f8fab10c72aec6e7a7d83762a40117a1e6e24','siAgeModel':'3e8e17567afa296213b414083eaca59b30455b81','PHARE':'0ab277b38ddac534d6d91e4580134cf6c78822ca'}
versions['PHARE_compatibility_environment']={'numpy':'1.26.4','scipy':'1.12.0','scikit-learn':'1.3.2','pandas':'2.2.3','celltypist':'1.6.3','torch':'2.2.2',
 'patch':'Optional autogenes import, unused by single-cell predict; original source preserved separately.',
 'status':'Model-only smoke, not full preprocessing validation.'}
(OUT/'versions_and_seed.json').write_text(json.dumps(versions,indent=2))
shutil.copy2(PROJECT/'outputs/model_manifest.json',OUT/'round1_model_manifest_preserved.json')

primary=[
 {'source':'HGNC current and historical archives','url':'https://www.genenames.org/download/archive/','role':'Approved stable HGNC identity, old names and Ensembl cross references'},
 {'source':'Ensembl BioMart current','url':'https://www.ensembl.org/biomart/martservice','role':'Successful GET response for residual model IDs; REST/archived failures retained as failures'},
 {'source':'GSE213516','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213516','role':'Independent candidate rejected: residual mask genes absent from feature universe'},
 {'source':'Author cell labels, GSE213516','url':'https://doi.org/10.6084/m9.figshare.23789844.v1','role':'Metadata readable; file download failures are not evidence of data absence or DUA'},
 {'source':'GSE96583','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583','role':'Eight control donor identities/labels verified; insufficient mask coverage, count download failed'},
 {'source':'Independent adult PBMC glaucoma atlas','url':'https://cellxgene.cziscience.com/collections/de2cde16-c8d3-4a6d-80be-1be9e879aaca','role':'Feature annotation inspected; no cell-level experiment performed'},
 {'source':'Independent 2025 immune aging atlas','url':'https://cellxgene.cziscience.com/collections/cc431242-35ea-41e1-a100-41e0dec2665b','role':'24-donor, age 20–75 study metadata; one ILC partition feature annotation inspected, not 24 matched-type validations'},
 {'source':'siAge reference','url':'https://www.synapse.org/Synapse:syn61609846','role':'Metadata and access-requirement endpoint readable; anonymous DOWNLOAD explicitly denied 403'},
 {'source':'PHARE author release','url':'https://doi.org/10.5281/zenodo.13838466','role':'H5AD HVG-only input rejected; archive index inspected, partial downloads not read as completed data'}]
(OUT/'evidence_sources.json').write_text(json.dumps(primary,indent=2,ensure_ascii=False))

table='\n'.join(f'| {r.celltype} | {r.missing_round1}/{r.nonzero_genes} | {r.missing_after_Ensembl_HGNC}/{r.nonzero_genes} | {100*r.residual_missing_fraction:.3f}% |' for r in g.itertuples())
report=f'''# 第2轮审计归档：按用户决定停留G1

状态：**PAUSED_AT_G1_BY_USER / 尚不能判定科学 PASS 或 KILL**。用户决定为“保持原门槛，停留 G1”。不开展探索性后续运行，也不自动恢复。固定种子 **{SEED}**。原始第1轮和 HIHA K0 归档均未改写。

## 步骤0：PASS

主分析冻结为 68 人、320 次稳态访视、90 个批次。按已固定配对规则，7–9 天、75–105 天、330–400 天分别有 33、28、62 人。320 次访视均已连接作者公开的整数年龄（25–67 岁），不由出生年倒算精确年龄。

与第1轮 89 人有 68 人重叠；与 87 人有 66 人重叠。K0 新纳入 BR2016、BR2047，因为两人在 K0 有4次其他合格稳态访视，而第1轮的候选 IV 配对均跨已记录接种。87 人中有21人不进入K0；逐人原因见 population_difference_donors.csv，其中5人疫苗日期不完整、BR2021访视标签歧义，其他人冻结90天窗口后不足4次合格访视。87人仅保留为敏感性定义。

## 步骤1：sc-ImmuAging覆盖率子项 PASS；G1 整体未完成

以下比例仅指sc-ImmuAging。siAge和PHARE尚未完成同等级的映射、完整预处理及遮蔽审计。

| 类型 | 第1轮缺失 | 当前残余缺失 | 当前残余比例 |
|---|---:|---:|---:|
{table}

映射依据包括稳定 Ensembl ID、HGNC 当前与 2022/2024 历史 HGNC 身份、HGNC 旧名/别名、NCBI org.Hs.eg.db 交叉引用及成功的 Ensembl BioMart GET 响应。歧义保留，不依据预测结果挑选映射。SOD2 的同一批准符号对应多个 GEO 特征时，映射表明确标记需先合并原始计数、再归一化；原始计数分母只能计算一次。并非已经完成320份矩阵的逐样本特征一致性核查。

八个固定 CELLxGENE 镜像的注释均已读取（32357个特征）；未发现需要因 >5% 残余缺失而切换主输入的情况。详细版本和字段见 CELLxGENE_annotation_inspection.json。

官方 CD8 示例：5名donor、1164个细胞，配对共享同一批100×15 pseudocell成员；连续 donor 均值预测的遮蔽相关 r={mask['Pearson_r']:.9f}，偏倚={mask['bias_years']:.6f}年，完整可用输入 MAE={mask['full_MAE']:.6f}年，遮蔽 MAE={mask['masked_MAE']:.6f}年，变化={mask['MAE_change']:.6f}年。按 r<0.9 规则，**CD8移出主分析**。示例自身只含592/631个非零权重特征，其余按官方函数填0；这不是独立队列准确度验证。

用户已明确同意：CD8 的可靠性、SEM及均值回归结果全部只放补充材料，不能支持最终KILL。90天主、7–9天敏感性的原选择仍记录，但CD8不再构成正式主判据。

非CD8独立验证输入尚未就绪。GSE213516的17人完整特征表不含残余遮蔽基因，不能把天然相同的完整/遮蔽预测作为G1通过。GSE96583虽含8名对照供体和作者类型标签，残余特征覆盖不足。两套新的CELLxGENE参考也只覆盖部分残余ID；G1_independent_validation_coverage.csv是**注释筛查**，不是已完成的匹配细胞类型遮蔽试验。未据此伪造 r、偏倚或 MAE。

## 步骤2–5：NOT_RUN

未越过尚未解决的G1开展68人主分析，故没有新的ICC、SEM、MDC95、LODO年龄加速、donor bootstrap CI、方差分解、RTM模拟或更新的新颖性结论。第1轮20人结果没有改名冒充68人结果。

siAge：匿名元数据可读、accessRequirement返回0条，但下载明确403：Anonymous users have only READ access permission。这不等于已证实需要DUA；未登录、未接受条款、未获取受控对象。完整原生标签转移尚未运行。

PHARE：已在兼容环境加载原始模型并对合成组成向量产生有限预测；这只属于工程 smoke test。未完成真实输入的全预处理运行验证，也未完成主队列可靠性。PHARE是供体层面细胞组成时钟，不能虚构为五个细胞类型的独立时钟。

## 步骤6：NOT_EVALUABLE

没有触发可报告的终局KILL，也不能宣布PASS。当前未知项是完整G1验证及下游必需证据；不能把未知写成可靠性低或均值回归成立。

## 可复现性

实际运行脚本包括 step0_population.py、export_orgdb_mapping.R、prepare_gene_mapping.py、read_mirror_annotation.py、read_sample_ages.py 和 G1_official_CD8_mask.R。build_checkpoint.py汇总状态。versions_and_seed.json、R sessionInfo、gene_annotation_receipts.json、K0_source_manifest.json和SHA256SUMS.csv记录版本、种子与哈希。实际运行过的脚本与仅下载的模型源代码严格区分。

此报告是未完成的检查点，不能作为完整第2轮审计结论、两钟验证或研究可发表性的证明。
'''
(OUT/'checkpoint_report_zh.md').write_text(report)

# Hash every user-facing file except the manifest itself. Source/input files are
# hashed separately without promoting partial archives to completed inputs.
manifest=[]
for p in sorted(OUT.rglob('*')):
 if p.is_file() and p.name!='SHA256SUMS.csv':
  manifest.append({'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(OUT/'SHA256SUMS.csv',index=False)
print('Checkpoint written; no unrun numerical results created.')
