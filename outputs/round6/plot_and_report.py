#!/usr/bin/env python3
"""Render frozen results only; no new models, correlations,or inferential tests."""
from pathlib import Path
import json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[2];O=P/'outputs/round6';CTS=['B','CD4T','CD8T','MONO','NK'];SEED=20261025
def md(x):
 x=x.copy().fillna('NA');return '| '+' | '.join(map(str,x.columns))+' |\n| '+' | '.join(['---']*len(x.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v).replace('|','/') for v in r)+' |' for r in x.itertuples(index=False,name=None))
def save(fig,name):
 for ext in ['png','svg','pdf']:fig.savefig(O/(name+'.'+ext),dpi=220,bbox_inches='tight')
 plt.close(fig)
def main():
 d=json.loads((O/'replication_decision.json').read_text());r=pd.read_csv(O/'replication_reliability_results.csv');cov=pd.read_csv(O/'model_gene_coverage.csv');pop=json.loads((O/'primary_population_freeze.json').read_text());kcounts=pd.read_csv(O/'K_c_qualifying_donor_counts.csv');interval=pd.read_csv(O/'external_primary_nominal_interval_distribution.csv');bs=pd.read_csv(O/'same_cross_batch_pair_counts.csv');allen=pd.read_csv(O/'Allen_published_vs_e_Pearson_r.csv');comp=pd.read_csv(O/'replication_comparison_with_AIDA.csv');qa=pd.read_csv(O/'all_type_count_reconstruction_QA.csv');aid=pd.read_csv(P/'outputs/round5/P4_external_technical_ICC.csv');aid=aid[aid.clock.eq('published')&aid.metric.eq('age_acceleration')].set_index('celltype');z=r[r.metric.eq('age_acceleration')].set_index('celltype').loc[CTS]
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False});fig,ax=plt.subplots(1,2,figsize=(11,4.5))
 for i,ct in enumerate(CTS):
  e=z.loc[ct];a=aid.loc[ct];ax[0].errorbar(e.ICC2_1,i-.1,xerr=[[e.ICC2_1-e.ICC_CI_low],[e.ICC_CI_high-e.ICC2_1]],fmt='o',capsize=3,color='#C57639',label=f'Terekhova longitudinal (n={pop["n_donors"]})' if i==0 else None);ax[0].errorbar(a.ICC2_1,i+.1,xerr=[[a.ICC2_1-a.ICC_CI_low],[a.ICC_CI_high-a.ICC2_1]],fmt='o',capsize=3,color='#7D8FA0',label='AIDA technical (n=20)' if i==0 else None);ax[1].errorbar(e.MDC95,i,xerr=[[e.MDC95-e.MDC95_CI_low],[e.MDC95_CI_high-e.MDC95]],fmt='o',capsize=3,color='#257976')
 for a in ax:a.set_yticks(range(5),CTS);a.invert_yaxis();a.grid(axis='x',alpha=.18)
 ax[0].axvline(.60,color='.6',ls='--',lw=.8);ax[0].set(xlabel='Published-clock age-acceleration ICC(2,1), 95% CI',title='A  Biological visits vs same-draw technical repeats');ax[0].legend(fontsize=7,frameon=False,loc='lower left');ax[1].set(xlabel='Age-acceleration MDC95 (years), donor bootstrap 95% CI',title='B  Longitudinal measurement resolution');fig.tight_layout();save(fig,'Figure6_external_longitudinal');z.reset_index().to_csv(O/'Figure6_external_plot_data.csv',index=False)
 fig,ax=plt.subplots(figsize=(7,3.8));ix=np.arange(5);ax.bar(ix-.18,allen.published_Pearson_r,.34,color='#7D8FA0',label='Published sc-ImmuAging');ax.bar(ix+.18,allen.e_Pearson_r,.34,color='#257976',label='Pseudobulk e (existing result)');ax.set(xticks=ix,xticklabels=CTS,ylabel='Donor-level Pearson r with chronological age',ylim=(0,1),title='Allen K0, 68 donors: existing CSV values only');ax.legend(fontsize=8,frameon=False);fig.tight_layout();save(fig,'Figure5_Allen_age_validity_readonly')
 short=[]
 for x in r.itertuples():short.append({'细胞':x.celltype,'指标':'预测年龄' if x.metric=='prediction' else 'LODO年龄加速','N':x.n_donors,'ICC [95%CI]':f'{x.ICC2_1:.3f} [{x.ICC_CI_low:.3f}, {x.ICC_CI_high:.3f}]','SEM [95%CI]':f'{x.SEM:.3f} [{x.SEM_CI_low:.3f}, {x.SEM_CI_high:.3f}]','MDC95 [95%CI]':f'{x.MDC95:.3f} [{x.MDC95_CI_low:.3f}, {x.MDC95_CI_high:.3f}]'})
 tab=md(pd.DataFrame(short));decision='PASS：复现' if d['replication_verdict']=='复现' else 'FAIL：未复现（完成有效分析）';steps=[{'步骤':'预设登记','状态':'PASS','证据':'2026-10-27T09:41:05.918983+00:00写入gate_ledger，先于下载/计算','一句话结论':'单一Terekhova来源、固定种子20261025、K关卡和复现阈值未改。'}, {'步骤':'K-a','状态':'PASS','证据':'五类HCA对象匿名GET200/206；完整文件SHA匹配；一般CCBY许可而非签署型DUA；已验证整数计数无损还原','一句话结论':'公开logCP10K加原始文库元数据可无损还原A1输入，不将log值当原始计数。'}, {'步骤':'K-b','状态':'PASS','证据':'166人/317生物样本单元；102人>=2访视；年度设计出处PDF第25页','一句话结论':'支持健康跨年度重复；本轮CSV未含实际日期，仅报告名义年度波次分布。'}, {'步骤':'K-c','状态':'PASS','证据':f'{pop["n_donors"]}人、{pop["n_visits"]}合格访视、{pop["n_cells"]:,}细胞；每类型每访视>=30','一句话结论':'主分析人群按元数据和nativeQC冻结，未按预测值筛选。'}, {'步骤':'预设外部复现','状态':decision,'证据':f'B残差ICC={d["B_age_acceleration_ICC"]:.6f}；MONO={d["MONO_age_acceleration_ICC"]:.6f}；低于AIDA技术重复{d["n_types_below_AIDA_published_technical"]}/5类','一句话结论':f'按两个AND条件的预设点估计规则判为{d["replication_verdict"]}；不以CI/p值替换阈值。'}, {'步骤':'Allen效度补表','状态':'PASS（只读既有结果）','证据':'读取round5/P2_chronological_age_signal.csv的已发表/e五类r','一句话结论':'未重新计算Pearson，未运行e/c。'}];pd.DataFrame(steps).to_csv(O/'step_status.csv',index=False)
 report=f'''# 单细胞免疫时钟可靠性审计：第6轮

本轮判定：**{d['replication_verdict']}**。固定种子 **20261025**，只运行已发表sc-ImmuAging。K关卡和复现判据于 **2026-10-27T09:41:05.918983+00:00** 写入台账，早于任何新下载/数据计算；分析人群于 **{pop['recorded_UTC']}** 冻结，早于任何预测。

{md(pd.DataFrame(steps))}

## K-a：匿名获取与真实矩阵格式

唯一数据源为 [HCA Terekhova健康血液图谱](https://explore.data.humancellatlas.org/projects/896f377c-8e88-463e-82b0-b2a5409d6fe4)。B包222,791,692字节匿名GET200、其余四类匿名范围读取GET200/206，之后完整对象SHA逐一匹配HCA公开记录。未使用登录、机构签字或申请型DUA。[HCA页面名为DataUseAgreement的通用说明](https://data.humancellatlas.org/about/data-use-agreement)仅列CC-BY4.0及release policy；这不等于没有许可条款。[AWS目录](https://registry.opendata.aws/singlecellhumanbloodatlas/)标注的Synapse controlled-access路线未用于本轮取数。

作者RNA输入是普通dgCMatrix log1pCP10K，而不是直接整数矩阵；但CSV保留逐细胞完整原始nCount_RNA及nFeature_RNA。对所有主分析细胞逐类验证 `expm1(x)*nCount_RNA/10000`：全部非零值恢复整数、列和匹配原文库总量、非零基因数匹配nFeature_RNA；最大误差见表。无原始总量或不能无损还原时不允许把归一化值当count。本轮没有这种替代。

{md(qa.round(12))}

CD4曾在HCA跳转端点出现连接超时和TLS EOF；改用HCA清单已列出的同一SHA256对象S3公开镜像，未换数据集，原失败receipt保留。作者代码SHA和官方文件对象记录均随归档保存。规范化数据在回归/Harmony之前载入，未用Harmony低维坐标代替表达矩阵。

## K-b / K-c：观察单位、时间和细胞覆盖

`Donor_id`为人，`Tube_id`为生物采样访视，`File_name`为技术文库；166人/317访视中64人一次、53人两次、49人三次，因此102人有重复采样。本轮取得的逐细胞CSV不含采样日期或相对天数；只提供整数年龄，不能把其差值当精确间隔。未穷尽所有原始临床文件。

年度采血设计在同一Terekhova数据的[原始再分析论文Methods](https://assets-eu.researchsquare.com/files/rs-6784948/v1_covered_9100f85d-e659-45c5-b3f5-61539c2feabd.pdf)PDF第25页明确描述；原HCA/AWS则描述三年收集期。访视波次按来源Tube编码解码：ABCDE为首轮、GHJKL第二轮、MNPST第三轮，F女性前缀单独剥除。**这是结合来源编码和年度设计的名义波次推断，不是原始采样日期字段。** 解码后每人每波次唯一、整数年龄单调；所有时间表保留actual_days为空。

应用A1常规谱系：B排除plasma；MONO仅classical/nonclassical，不纳入DC；CD8排除含混NKT-like；CD4/NK保持来源常规大类。保留200<=nFeature_RNA<=5000、percent.mt<=25、原始文库>0的细胞。五类交集313次访视满足每类>=30；其中至少有两次合格访视的101名donor共250次访视、1,421,481个细胞进入分析。

{md(kcounts)}

每人选最早与最晚合格波次构成101对；LODO使用这101人的全部250次共同合格访视，含可用的中间访视，各donor总权重相等。单访视供者不加入年龄校准。名义间隔分布：

{md(interval)}

因此这是约一年/约两年混合的纵向重测，不能表述为精确365天或730天分析。逐类批次配对数为：

{md(bs)}

97对跨Batch、4对同Batch在五类一致；多批次共享或未知类别为0。File_name共享情况另存台账，不能把它当新的生物访视。

## 基因映射和复现结果

RNA特征为gene symbol，经冻结的Ensembl模型ID→HGNC批准名/历史别名及模型原始注释链接到来源特征；多义别名不择第一行。同一批准基因多来源特征先在CP10K/count尺度求和、再log1p。未测/未解析特征依native补零并逐基因列权重；没有拿工程缓存缺失当生物零。

{md(cov[['celltype','nonzero_model_genes','mapped_genes','residual_missing_genes','residual_missing_fraction']].round(6))}

原始权重和截距冻结，不训练或调参。对每donor×visit×类型从合并的技术文库细胞中生成100个15细胞pseudocell，使用A1连续均值、不取整。原始预测年龄及等供者权重LODO年龄加速均保留。2000次共享donor聚类bootstrap采用SeedSequence[20261025,6201]，每次重拟合LODO并排除被留出原始donor的全部副本；直接从每个指标分布取百分位95%CI。SEM使用配对观测SD×sqrt(1−ICC)，MDC95=1.96×sqrt(2)×SEM，单位年。所有bootstrap有效数和数值参考检查均保存。

{tab}

与已发表模型AIDA同次采血技术重复的比较：

{md(comp[['celltype','ICC2_1','AIDA_technical_ICC_published','strictly_below_AIDA_technical']].round(6))}

预设第一条腿B和MONO残差ICC均<=0.60：**{d['B_and_MONO_both_le0_60']}**；第二条腿至少3/5类低于AIDA技术ICC：**{d['at_least3_types_below_AIDA']}**，实际 **{d['n_types_below_AIDA_published_technical']}/5**。两个AND条件最终判为 **{d['replication_verdict']}**。未改变阈值、改选访视、换数据集或使用e/c来救回结论。

该判定针对已报告注释映射/native缺失策略下的时钟应用。跨队列技术、细胞数、时间和人群有差异；97/101对跨批次，故不能将可靠性差异因果归为纯生物衰老。AIDA仅技术重复，Terekhova为年度健康访视，两者estimand不同；低于技术重测是预设操作性复现条件，并非独立的因果证明。

## Allen图5效度补表：只读已有CSV

{md(allen[['celltype','n_donors','published_Pearson_r','e_Pearson_r']].round(6))}

读取round5已有Pearson结果，未重算、未运行e/c。MONO两模型r均较低，NK的e也没有高于已发表模型。此效度信息与重测ICC分别解释。

## 审计与预注册状态

预注册P2的验证性/探索性标签继续 **待声明**；既有CD4基因ICC/MSR/MSC/MSE接触声明保留。本轮没有新增CD4基因层面方差成分分析，只做输入完整性及预设时钟可靠性。

脚本、种子、版本、输入/输出SHA256和checkpoint见本目录。图6展示外部残差ICC/MDC及AIDA参考；图5新增panel仅画Allen已有r。大矩阵留在可恢复工作缓存，checkpoint不重复打包数GB来源数据。
'''
 (O/'round6_report_zh.md').write_text(report)
 results='；'.join(f'{ct}：{z.loc[ct,"ICC2_1"]:.3f}（95%CI {z.loc[ct,"ICC_CI_low"]:.3f}–{z.loc[ct,"ICC_CI_high"]:.3f}），SEM {z.loc[ct,"SEM"]:.3f}年，MDC95 {z.loc[ct,"MDC95"]:.3f}年' for ct in CTS)
 methods=f'''# 第6轮：方法与结果更新段落

方法：在下载和计算前冻结匿名获取、至少15名健康重复供者、每类每访视至少30个细胞的K关卡，以及B和MONO年龄加速ICC点估计均≤0.60且至少3/5类低于AIDA已发表模型技术重复ICC的联合复现标准。仅使用Terekhova健康血液图谱的HCA公开对象；归一化矩阵配合原始文库总量经逐值整数及列和验证无损还原计数，使用已发表sc-ImmuAging固定权重、HGNC历史映射、native缺失策略与100×15细胞pseudocell连续均值。101名供者的250次合格访视进入LODO年龄校准；每供者以最早和最晚合格波次形成一对，共101对。名义年度波次间隔为一年39对、两年62对；本轮CSV未含实际采样日期，实际日间隔未取得；Tube编码波次为来源设计支持的推断，不将整数年龄差替代实际间隔。五类均有97对跨批次、4对同批次。采用2000次供者聚类bootstrap并每次重新拟合LODO，报告ICC(2,1)、SEM和MDC95及百分位95%CI。未运行变体e/c或训练新模型。

结果：原始年龄和LODO年龄加速完整指标列于复现结果表。年龄加速结果为{results}。B和MONO均≤0.60条件为{d['B_and_MONO_both_le0_60']}，低于AIDA技术重复为{d['n_types_below_AIDA_published_technical']}/5类，按冻结判据判为“{d['replication_verdict']}”。这一操作性结果反映年度健康访视条件下的重测表现；它受到来源特征映射、细胞数、跨批次和时间结构共同影响，不能据此将访视差异全部归为生物衰老。

效度补充：Allen68名供者的既有CSV显示，已发表/e的年龄相关r分别为B0.086/0.605、CD40.813/0.888、CD80.820/0.828、MONO0.158/0.169、NK0.344/0.325。本轮只读取这些结果用于图5，不作新的精度计算。预注册P2标签仍为待声明，先前CD4基因方差信息接触需如实披露。
''';(O/'manuscript_methods_results_zh.md').write_text(methods)
if __name__=='__main__':main()
