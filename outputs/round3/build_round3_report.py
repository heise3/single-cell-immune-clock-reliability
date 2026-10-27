#!/usr/bin/env python3
"""Assemble manuscript-facing tables and a bounded interpretation of frozen F."""
from pathlib import Path
import json
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3';A=P/'outputs/round2_A1'
def link(name,label=None):return f'[{label or name}]({O/name})'
def md(df):
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in df.itertuples(index=False,name=None))+'\n'
def main():
 f=json.loads((O/'R2_F_decision.json').read_text());r=pd.read_csv(O/'R2_E_reliability_results.csv');a=pd.read_csv(O/'R2_accuracy_results.csv');g=pd.read_csv(O/'R2_F_gate_by_type.csv');models=pd.read_csv(O/'R2_model_training_results.csv');old=pd.read_csv(A/'reliability_results.csv');rt=pd.read_csv(O/'R1_critical_SEM_vs_CD8_estimates.csv');lit=pd.read_csv(O/'R3_overlapping_literature_2025_2026.csv')
 w=r.pivot(index=['celltype','stratum','n_paired_donors'],columns='variant',values=['ICC2_1','ICC_CI_low','ICC_CI_high','delta_ICC_vs_replicated_baseline','delta_ICC_CI_low','delta_ICC_CI_high']);w.columns=[f'{v}__{stat}' for stat,v in w.columns];w.reset_index().to_csv(O/'R2_manuscript_reliability_wide.csv',index=False)
 w=a.pivot(index=['celltype','variant'],columns='dataset',values=['n_donors','MAE','Pearson_r']);w.columns=[f'{ds}__{stat}' for stat,ds in w.columns];w.reset_index().to_csv(O/'R2_accuracy_wide.csv',index=False)
 gi=pd.read_csv(O/'R2_S_gene_ICC.csv');gs=gi.groupby('celltype').agg(n_genes=('feature_id','size'),n_estimable=('ICC2_1','count'));gs['ICC_ge_0_5']=gi.groupby('celltype').ICC2_1.apply(lambda z:int(z.ge(.5).sum()));gs['ICC_ge_0_7']=gi.groupby('celltype').ICC2_1.apply(lambda z:int(z.ge(.7).sum()));gs.reset_index().to_csv(O/'R2_gene_ICC_filter_summary.csv',index=False)
 models['lambda_at_frozen_minimum']=np.isclose(models.lambda_fraction,1e-4);models['baseline_K_at_frozen_2000_cap']=models.baseline_selected_K.eq(2000);models[['celltype','variant','baseline_selected_K','lambda_fraction','lambda_at_frozen_minimum','baseline_K_at_frozen_2000_cap']].to_csv(O/'R2_model_grid_boundary_flags.csv',index=False)
 scope=[]
 for ct in ['B','CD4T','CD8T','MONO','NK']:
  z=old[old.clock.eq('sc-ImmuAging')&old.celltype.eq(ct)&old.metric.eq('age_acceleration')&old.stratum.eq('75-105d')].iloc[0];z2=r[r.celltype.eq(ct)&r.variant.eq('a')&r.stratum.eq('75-105d')].iloc[0];scope.append({'celltype':ct,'A1_published_model_K0_68_90d_ICC':z.ICC2_1,'A1_paired_donors':z.n_donors,'R2_retrained_OneK_baseline_E_90d_ICC':z2.ICC2_1,'R2_paired_donors':z2.n_paired_donors,'interpretation':'Different models and evaluation populations; difference is NOT a gene-repair effect and is never used by F'})
 pd.DataFrame(scope).to_csv(O/'R0_vs_R2_model_and_population_scope.csv',index=False)
 step=pd.DataFrame([{'step':'F pre-registration','status':'PASS','evidence':'gate_ledger.json round3_F recorded2026-10-27T00:54:41.574244+00:00 before calculation; F_prespecified.json retained','one_sentence_conclusion':'F门槛和比较对象未因结果更改。'},{'step':'R0','status':'PASS','evidence':'15 wide rows;42 variance rows;batch largest31/42;all5 scImmu adjusted rows donor largest','one_sentence_conclusion':'汇总完成；整体批次最多，sc-ImmuAging五类内donor均为最大单项。'},{'step':'R1','status':'PASS','evidence':'56 SEM grid points;100000 simulations each;Holm critical SEM1.8y,unadjusted2.0y','one_sentence_conclusion':'90天与一年SEM的CI下限高于临界点；短间隔CI跨临界点。'},{'step':'R2','status':f['repair_verdict'],'evidence':f"25 models;S34/E34;S short pairs14;E90d13/year31;2000 valid donor bootstrap each;qualifying variants={','.join(f['successful_variants'])}",'one_sentence_conclusion':'c和e满足F；修复结果可写入稿件，但准确性与外推限制必须保留。'},{'step':'R3','status':'PASS','evidence':'OpenAlex5 successful;S2 ordinary3 queries429 retained;bulk4 queries successful;13 curated partial overlaps;no full overlap identified in inspected evidence','one_sentence_conclusion':'完成有界API补查；未发现完全重叠，不等于穷尽排除。'}]);step.to_csv(O/'step_status.csv',index=False)
 tt=g[g.variant.isin(['c','e'])&g.type_meets_F].copy();tt['baseline90']=tt.ICC2_1-tt.delta_ICC_vs_replicated_baseline;summary=pd.DataFrame({'变体':tt.variant,'类型':tt.celltype,'基线ICC':tt.baseline90.map(lambda x:f'{x:.3f}'),'变体ICC [95%CI]':[f'{v:.3f} [{lo:.3f}, {hi:.3f}]' for v,lo,hi in zip(tt.ICC2_1,tt.ICC_CI_low,tt.ICC_CI_high)],'ΔICC [95%CI]':[f'{v:+.3f} [{lo:+.3f}, {hi:+.3f}]' for v,lo,hi in zip(tt.delta_ICC_vs_replicated_baseline,tt.delta_ICC_CI_low,tt.delta_ICC_CI_high)],'OneK ΔMAE/岁':tt.MAE_increase.map(lambda x:f'{x:+.3f}')})
 report=f'''# 单细胞免疫时钟可靠性审计：第3轮成稿分析

**预设修复门槛 F：{f['repair_verdict']}。** 变体 c（基因 ICC≥0.7）和 e（donor pseudobulk）均在 MONO、NK 两类达到约90天年龄加速 ICC 提升≥0.15，且相应 OneK1K 留出集 MAE 增加≤1.5年。可采用“可靠性审计＋预设修复比较”的稿件结构；这不等于模型已具备准确、可迁移的绝对年龄标定。固定种子 **20261025**。

{md(step[['step','status','one_sentence_conclusion']])}

## 预设、恢复与证据边界

F 于 **2026-10-27 00:54:41 UTC**（本地2026-10-26）先写入活动 gate_ledger.json，再开展本轮计算；{link('F_prespecified.json')}保留原文。基准是同等 OneK1K 训练条件下的复刻模型 a，F 的两项要求须由**同一变体、同一细胞类型**同时满足，且至少两类属于 B/MONO/NK。未采用旧发表模型或完整68人结果作为修复对照。

恢复时核验旧 round2_A1_checkpoint.zip 的 SHA256 为 `9342f9f9562dc9c70d71daf1755e626d6f081fb05a13f6373c5dc06bdbc1ff7e`，并逐项核对其272个清单文件；{link('checkpoint_recovery_QA.json')}。旧台账、原结果和归档保持不变。R0只重排既有结果，没有重算可靠性和方差分解。

## R0：PASS

{link('Table1_scImmuAging_reliability_wide.csv','表1')}含五类×三个间隔的15行，每行同时列出原始预测与LODO年龄加速的 ICC、SEM、MDC95及95%CI。主队列68人、320次稳态访视；三个间隔配对数依次为33、28、62。{link('Table2_variance_components_42.csv','表2')}包含42项各方差成分占比及最大项。

**一句话结论：合并42项时，批次在31项中是最大单项；sc-ImmuAging的原始/年龄调整评分内，donor在五类中均为最大单项。** RAM Up原始评分有1类donor最大，其余RAM行均为批次最大；不能把总体31/42误写成sc-ImmuAging的批次主导。最大单项也不必占总方差50%以上。

“访视”列实际为**访视嵌套＋未建模残差**。一个访视只有一个最终评分，无法再区分自然生物波动与未建模处理误差。细胞抽样和pseudocell方差来自A1固定的各100次条件重采样；占比不是因果归因。RAM未标定为岁，不报MDC年数。

## R1：PASS

38名配对CD8变化沿用Fig.5源数据：全队列均值+1.222845岁，负变化16/38，负变化亚组中位数−2.753863岁。SEM从0.5至6.0岁、步长0.1，共56点，每点100,000次模拟，零假设为 Δ~N(全队列观测均值, √2·SEM)。两项双侧经验p按第2轮规则计算，并对负变化人数、负子组中位数做Holm(2)校正。

**任一统计量校正p<0.05的最大网格SEM为1.8岁；未校正敏感性为2.0岁。** 这是网格临界值，不插值，不是全56点的同时置信带。模拟使用同一组标准正态随机数减小曲线蒙特卡洛抖动。负子组为空时仅从中位数统计量中排除，相应有效模拟数已记录。

{md(pd.DataFrame({'间隔':['7–9天','75–105天','330–400天'],'CD8残差SEM [95%CI]/岁':['2.081 [1.620, 2.503]','3.377 [2.382, 4.265]','3.536 [2.671, 4.439]'],'SEM/临界SEM':['1.156','1.876','1.965'],'CI下限/临界SEM':['0.900','1.323','1.484']}))}

所有原始/残差口径、三个间隔及CI上下限的双向比值见{link('R1_critical_SEM_vs_CD8_estimates.csv')}；完整曲线见{link('R1_SEM_sensitivity_curve.csv')}。**一句话结论：90天与一年SEM的CI下限仍高于1.8岁，而7–9天CI下限落在临界点下方。** 指定零假设以观测全队列均值为中心，检验所选负向人数/幅度是否异常；不检验平均BCG效应为零，未拒绝也不能证明均值回归机制或BCG无效。将HIHA SEM用于BCG仍依赖跨研究可迁移假设。

![R1临界点]({O/'Figure_R1_SEM_tipping_point.png'})

## R2：PASS（预设F）

K0按既有年龄带（<45/≥45岁）分层随机对半，S=34、E=34；S含14名7–9天配对donor。E含165次访视，90天配对13人，一年配对31人。S只提供基因ICC，不参与时钟训练；E从未用于本轮特征选择、超参数或PCA维数选择。**E已在A1描述性审计中出现，本轮拆分提供训练泄漏防护，不构成全新的、未接触过的确认队列。** {link('R2_K0_donor_split.csv')}、{link('R2_E_pair_ledger.csv')}。

OneK1K固定Census2025-11-08和源版本，981名donor，784训练/197留出；MONO因源数据覆盖为783/195。五类计数在native QC后为1,204,490细胞。健康AIDA固定v2排除6名LONZA商业共同对照，剩619名、19–77岁，五类native QC后1,018,634细胞；并非每人都有五类，逐类型N均单独报告。全球is_primary_data标记的v1/v2重复未用于过滤单独选定v2内的细胞；没有合并v1。三个计数来源的每个完整raw行和均与Census raw_sum完全一致。

训练候选为OneK与HIHA共同注释的31,750个Ensembl基因，不使用发表模型的年龄筛选特征清单。S的donor×visit×类型pooled pseudobulk取log2(1+CPM)，仅用7–9天配对计算基因ICC。a为训练折内Pearson年龄筛选、BH q<0.05、100–2000基因数网格、LASSO；b/c与S基因ICC≥0.5/0.7取交集；d训练折内PCA后LASSO；e donor pseudobulk输入。完整方法预设见{link('R2_execution_contract.json')}。

**复刻范围须明示：** sc-ImmuAging最终模型是弹性网的α=1端点LASSO。本轮复刻OneK-only Pearson筛选分支，并加入donor阻断CV与折内筛选；没有重建作者五队列、MIRP/MIRA等全部训练过程。a–d保留log1p CP10K与100×15 pseudocell、逆年龄频数权重、无基因标准化。donor使用连续均值，与A1一致；发表正文中位数及包内均值取整均非本轮精度口径。e同时改变聚合方式与表达尺度至log2(1+CPM)，不能将其效应仅归于取消pseudocell。{link('R2_protocol_replication_QA.csv')}。

五折donor CV选择a的基因数与λ；b–e保持同一基因数规则并重调λ。所有变体用同一OneK留出donor和E配对，2000次E donor聚类bootstrap，每次重新拟合LODO。全部50项ICC均有2000次有效bootstrap。CIs条件于已经训练完成的固定模型，不含重训的不确定性。负ICC不截零。

{md(summary)}

c、e各有2类满足F；b为0类，d为1类。**F依据点估计预设判定，未要求ΔICC的CI下限≥0.15；MONO两项ΔICC的CI跨零，不能改写成已显著确认提升。** NK的点提升接近0.15门槛，也应保留数值精度。B的90天ICC没有被c/e改善。所有五类与一年结果均见{link('R2_manuscript_reliability_wide.csv')}，不只报告成功类型。

![R2可靠性]({O/'Figure_R2_E_reliability.png'})

**外部精度限制：** c在AIDA的MONO/NK MAE分别25.99/24.92岁，Pearson r为0.162/0.094；e分别14.72/21.95岁，r为0.288/0.332。a的AIDA B MAE达到55.05岁。新模型非零基因在E与AIDA源注释覆盖均100%，因此不能将这些误差归为注释缺基因；跨队列分布偏移、建模及绝对标定问题仍需面对。本轮不对E/AIDA作结果驱动的平移、校准或重训。{link('R2_accuracy_results.csv')}提供75项MAE/r及各自donor数，{link('R2_prediction_calibration_diagnostics.csv')}提供偏倚/预测范围。

OneK留出20% c的MONO/NK MAE为14.77/13.57岁，e为15.35/12.39岁，满足的只是**相对a**误差增幅门槛。高ICC不能代替年龄准确性、临床有效性或外部纵向验证。{link('R0_vs_R2_model_and_population_scope.csv')}单列A1发表模型与本轮OneK复刻的口径差异，二者不作修复前后对照。

优化过程没有扩展结果边界：部分模型λ在预设下限、CD4基因数在2000上限，均原样报告{link('R2_model_grid_boundary_flags.csv')}。LASSO坐标下降耗时后，以等价LARS路径加速CD4/CD8，仍在原100个λ上取解，并逐拟合核查对偶间隙；真实训练样本与原生glmnet最大预测差约0.000053年，2000基因CD8训练折对偶间隙约1.55×10⁻¹⁰。{link('R2_solver_engineering_update.json')}保留改变和事前脚本哈希，未改变F。

**一句话结论：修复比较按预设F通过，证据适合支持有边界的可靠性改进，不适合宣称得到通用的高精度免疫年龄时钟。**

## R3：PASS（有界检索）

OpenAlex五条查询均成功；普通Semantic Scholar中三条持续429，错误回执保留；随后四组Boolean bulk查询均200，命中数为3/0/0/0，三组阳性对照7/2/29命中。普通结果548条＋bulk3条，标题去重516；有界首批查询每条最多100项，不是穷尽系统综述。逐条候选与读取层级见{link('R3_overlapping_literature_2025_2026.csv')}，API回执和筛选优先级表一并保存。

最强新增部分重叠为[Mishra等2026](https://www.nature.com/articles/s41467-026-73218-6)的333人纵向bulk血液基因ICC；[Sehgal等2026](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.70635)的DNAm生物/技术可靠性；[Borrus等2026](https://onlinelibrary.wiley.com/doi/10.1111/acel.70756)的干预假阳性；[TranslAGE平台](https://www.translage.io/)及[2026干预响应研究](https://www.nature.com/articles/s41591-026-04562-9)。[Park等2026](https://www.nature.com/articles/s41467-026-76737-4)使用重叠PBMC队列做性别分层年龄建模，亦属部分重叠。

13项部分重叠记录中，2项仅核查到主来源索引摘要，全文重叠仍未核实；普通S2访问限制及首批排名边界不能被“0命中”抹掉。**一句话结论：已审证据未发现完全重叠；可主张的增量限定为单细胞免疫donor重测、MDC/噪声分解、选择后BCG统计与匹配修复比较，可靠性审计、基因ICC筛选或PC去噪本身已有先例。**

## 交付与复现

表1、表2、R1曲线、R2全部ICC/精度/门槛表、25个新模型、逐基因权重、脚本及执行日志均在本目录。{link('REPRODUCE.md')}列运行顺序；{link('versions_and_seed.json')}固定版本与随机子流；{link('SHA256SUMS.csv')}覆盖交付文件；input_files_SHA256.csv及原始下载manifest绑定大型输入。原始计数、设计矩阵和版权全文缓存不重复打包。恢复与重跑需要原A1 checkpoint及README指定的公开缓存/源版本。正文可使用{link('manuscript_results_methods_zh.md')}中的方法与结果段落，须保留本报告中的限制。
'''
 (O/'analysis_report_zh.md').write_text(report)
 methods='''# 可用于成稿的方法与结果段落

方法：我们基于预先归档的健康人稳态访视台账开展供者层面的单细胞免疫年龄时钟可靠性审计。主队列为68名供者、320次访视；在7–9、75–105及330–400天窗口内，每名供者按预设距离和平局规则选取一对访视。原始预测与leave-one-donor-out年龄残差分别计算双向随机、绝对一致性单次测量ICC(2,1)。年龄残差模型含年龄及截距，供者总权重相等，排除被评估者的所有访视。SEM为配对测量合并SD乘以√(1−ICC)，MDC95为1.96√2·SEM；CI由2000次供者聚类bootstrap给出，残差模型在各次bootstrap中重新拟合。访视方差与未建模误差无法分离；细胞抽样与pseudocell随机性采用条件重采样方差。

BCG敏感性：沿用已公布源数据中38名供者的CD8配对变化，模拟以观测全队列均值为中心的正态变化，标准差为√2·SEM。SEM按0.1年步长从0.5至6.0年，每点100000次，比较负变化人数与负变化子组中位数，并作双侧经验检验和Holm(2)校正。该检验针对选择后统计量，不检验全队列平均效应为零。所有SEM点使用共同随机数，临界点仅按离散网格定义。

修复比较：68名供者按年龄带分层随机分为34名筛选供者和34名评估供者。筛选组14人的短间隔配对用于donor×visit×细胞类型pseudobulk log2(1+CPM)基因ICC；时钟训练仅使用OneK1K。OneK1K按供者分出20%固定留出集，其余五折供者阻断CV内重新进行Pearson年龄筛选、BH校正和参数选择。共同注释基因全集不依赖年龄效果；未使用发表模型的年龄特征清单。复刻基线为OneK-only Pearson筛选分支的LASSO，采用100个15细胞pseudocell、log1p CP10K、逆年龄频数权重和无基因标准化，并以连续供者均值为评分。其余变体为短间隔基因ICC≥0.5、≥0.7、训练内PCA和donor pseudobulk。此复刻不等于作者五队列完整模型；donor pseudobulk变体同时改变表达尺度。

评估：评估组有13人的90天配对与31人的一年配对。每个固定模型的LODO残差ICC采用2000次供者聚类bootstrap；比较使用同一抽样以估计配对ΔICC的CI。OneK1K留出集、Allen评估组及619名健康AIDA供者用于MAE和Pearson相关；未据其表现选择模型或校准。预设F要求同一变体在B/MONO/NK至少两类同时满足90天残差ICC增幅≥0.15及相应OneK留出MAE增幅≤1.5年。门槛按点估计判定，未预设显著性要求。评估供者曾用于此前描述性审计，故本拆分是本轮防训练泄漏措施，不是新独立确认研究。

结果：42项原模型方差分解中31项以批次为最大单项，sc-ImmuAging五类原始及年龄调整评分均以供者为最大单项。BCG模拟的Holm校正临界SEM为1.8年；90天CD8残差SEM为3.377年（95%CI 2.382–4.265），一年为3.536年（2.671–4.439），短间隔为2.081年（1.620–2.503）。仅短间隔CI覆盖临界点以下范围。此结果表明指定负向选择统计与所估重访波动可以相容，并不证明BCG平均效应为零或均值回归为唯一机制。

修复结果：基因ICC≥0.7变体使MONO的90天残差ICC由0.408提升至0.696、NK由0.804至0.961；OneK留出MAE分别增加0.525及1.216年。Donor pseudobulk使相同两类ICC达到0.761和0.957，MAE分别增加1.101及0.041年。二者达到预设F；其余变体未在两类同时达到。MONO提升的bootstrap CI跨零，90天样本仅13人，CI条件于固定训练模型，因此结果应视为初步修复证据。AIDA上相应MONO/NK模型精度仍有限，不能由高ICC推断正确绝对年龄或临床意义。

讨论定位：既有DNAm时钟可靠性、PC去噪、干预假阳性和纵向血液基因ICC研究构成明确先例。本研究的贡献应限定为公开单细胞免疫时钟在冻结稳态供者/窗口下的重测误差、年龄残差与噪声分解，以及匹配对照的修复验证。API检索有界且部分候选全文未核实，不作穷尽性首创声明。
'''
 (O/'manuscript_results_methods_zh.md').write_text(methods)
 readme='''# Round3 reproduction

Seed:20261025. Run from the project root with `work/analysis_env/bin/python`. Retain the immutable `outputs/round2_A1/` tree restored from the SHA-bound A1 checkpoint; new scripts import its LODO/ICC routines and metadata readers. Paths of original sources remain in the manifests. Do not treat metadata sparse range-read files as expression matrices.

Order actually executed:

1. Record F in `outputs/gate_ledger.json` before calculations; retain `F_prespecified.json`.
2. `run_R0_R1.py`: verify A1 archive/files, reshape R0 and simulate R1.
3. `preflight_Census.py`, `freeze_R2_contract.py`: inspect fixed source versions, freeze contract/donor splits before gene/model outcomes.
4. `aggregate_training_inputs.py`: full-count libraries/QC; OneK1K and S aggregates; S gene ICC. `prepare_training_design.py`: fold-local training ranks, exact native100×15 pseudocell memberships and designs.
5. `read_AIDA_source_metadata.py` uses `metadata_range_reader.py` to recover only needed source metadata. `prepare_evaluation_inputs.py` builds E/AIDA inputs using full raw library denominators. This input preparation did not predict E outcomes.
6. `train_repair_models.py`: native validation via `validate_Gram_native.R`, then25 models. Completed B/MONO/NK use exact Gram coordinate descent; CD4/CD8 use equivalent LARS with dual-gap checks and the same lambda grid. Known partial coordinate runs were interrupted, not accepted. Final script supports either all types or environment `R3_TRAIN_TYPES=CD8T`, `R3_BLAS_THREADS=2`. Completed-type checkpoints are reused, so move only `work/round3/models/` to a separate historical cache if intentionally performing a complete fresh rerun. Do not delete published artifacts. Deterministic rerun may differ by floating-point tolerances across solver/runtime versions.
7. `evaluate_repair_models.py`: no training; fixed-model accuracy, E LODO and2000 shared donor bootstrap draws, unchanged F.
8. `search_OpenAlex_SemanticScholar.py`, `refine_SemanticScholar_bulk.py`, `validate_S2_search_controls.py`, `finalize_R3.py`: bounded APIs, preserved429s and curated source review. Primary body review is documented, not automatically recreated by title matching.
9. `plot_round3.py`, `build_round3_report.py`, `package_round3.py`: figures, report, independent artifact checks, old-input hashes and archive.

Large caches in `work/round3/` are not in the deliverable ZIP. Raw count manifests retain fixed release, chunk hashes and size receipts; input_files_SHA256.csv distinguishes freshly rehashed artifacts from download hashes with size rechecks. New25 model NPZs and effective Ensembl coefficients are bundled. The NPZ prediction formula is `prediction = intercept + sum(input * coefficient)`; its expression_unit must be honored. All input normalization uses full released raw library sums, never only selected genes. PCA weights are folded into effective gene coefficients, so original PCs are unnecessary to predict the frozen linear models.

R0 CSVs are exact reshapes of old completed scientific results. R1 is the only repeated BCG null simulation. R2 is a distinct matched OneK-only model comparison, not a rerun of the published pretrained model. F is point-estimate based; wide CIs and poor external absolute accuracy must remain in any manuscript.
'''
 (O/'REPRODUCE.md').write_text(readme);print('REPORT_COMPLETE',f['repair_verdict'])
if __name__=='__main__':main()
