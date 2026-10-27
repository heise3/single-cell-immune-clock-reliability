#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';R=P/'outputs/round3'
def link(f,label=None):return f'[{label or f}]({O/f})'
def md(x):return '| '+' | '.join(x.columns)+' |\n| '+' | '.join(['---']*len(x.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in x.itertuples(index=False,name=None))+'\n'
def main():
 dec=json.loads((O/'F2_decision.json').read_text());d=pd.read_csv(O/'Figure4_F2_paired_delta_ICC.csv');t=pd.read_csv(O/'T_REML_batch_components.csv');tc=pd.read_csv(O/'T_batch_contribution_by_type_interval.csv');c=pd.read_csv(O/'Figure4_c_repeated_crossfit_summary.csv');a=pd.read_csv(O/'Figure5_accuracy_bias_decomposition.csv');bias=pd.read_csv(O/'AIDA_equal_type_error_summary.csv');cov=pd.read_csv(O/'AIDA_a_c_e_feature_coverage.csv')
 wide=a.pivot(index=['celltype','variant'],columns='dataset',values=['n_donors','MAE','Pearson_r','bias_mean_error']);wide.columns=[f'{ds}__{stat}' for stat,ds in wide.columns];wide.reset_index().to_csv(O/'Figure5_accuracy_wide.csv',index=False)
 status=pd.DataFrame([{'step':'F2预设','status':'PASS','evidence':'2026-10-27T02:28:44.808874+00:00先写gate_ledger；e主要、c次要','conclusion':'主要变体与F2未随结果更改。'},{'step':'T1来源拆解','status':'PASS','evidence':'31项均RAM：Down原始8/调整8；Up原始7/调整8；合并320访视、非间隔模型','conclusion':'不能把31/42写成五类sc-ImmuAging或某个特定间隔的批次主导。'},{'step':'T2同/跨批次','status':'PASS','evidence':'7–9天18/15对；90天5/23对；一年0/62对；n<3不估计','conclusion':'分层完成，但一年没有同批次参照，90天同批次仅5人。'},{'step':'T3批次调整','status':'PASS','evidence':'五类×三间隔ICC/SEM/MDC及条件bootstrapCI；其他donor批次BLUP，保留访视残差','conclusion':'完成允许的评分层面混合模型去批次；没有运行ComBat-seq或宣称RNA矩阵已整合。'},{'step':'T4贡献解释','status':'PASS','evidence':'批次关联占非donor方差32.7–64.8%；严格技术/生物分离不可识别','conclusion':'给出模型依赖的近似贡献，不把剩余波动全部称为生物学。'},{'step':'F2-e完整队列','status':dec['status'],'evidence':'68人，28/62配对；2000次配对donorbootstrap，4项目标差值CI>0','conclusion':'按预设F2内部确认，MONO和NK两个间隔均确认；B未确认。'},{'step':'F2-c交叉拟合','status':'PASS','evidence':'50×2折×5类=500模型；每训练半队列12–21短间隔配对；全50次均可估计','conclusion':'次要结果完成；50次百分位是拆分波动范围，不是bootstrapCI。'},{'step':'精度与AIDA偏倚','status':'PASS','evidence':'a/c/e三队列45项MAE/r/偏倚；AIDA五类等权偏移平方份额a86.6%,c62.4%,e47.7%','conclusion':'误差来源随模型/类型变化，移除均值偏倚不会改变Pearson相关。'},{'step':'B采样间隔','status':'PASS','evidence':'两篇原始Methods均T0至T3m；匹配75–105天SEM3.3768，临界比1.876/CI下限比1.323','conclusion':'三个月层是时间尺度匹配；个体准确日期未释放，流程可迁移仍需假设。'}]);status.to_csv(O/'step_status.csv',index=False)
 target=d[d.F2_target_type].copy();tab=pd.DataFrame({'类型':target.celltype,'间隔/天':target.stratum,'配对数':target.n_pairs,'a ICC':target.baseline_a_ICC.map(lambda z:f'{z:.3f}'),'e ICC':target.variant_e_ICC.map(lambda z:f'{z:.3f}'),'ΔICC [配对bootstrap95%CI]':[f'{v:+.3f} [{lo:+.3f}, {hi:+.3f}]' for v,lo,hi in zip(target.delta_ICC_e_minus_a,target.delta_CI_low,target.delta_CI_high)],'F2':target.CI_confirms_improvement.map({True:'确认改善',False:'未确认'})})
 tt=t.merge(tc[tc.stratum.eq('75-105d')],on='celltype',validate='one_to_one');techtab=pd.DataFrame({'类型':tt.celltype,'批次/总方差':tt.batch_fraction_total.map(lambda z:f'{100*z:.1f}%'),'批次/非donor方差':tt.batch_fraction_non_donor_variance.map(lambda z:f'{100*z:.1f}%'),'90天ICC调整前→后':[f'{b:.3f} → {v:.3f}' for b,v in zip(tt.ICC_before,tt.ICC_batch_adjusted)],'调整后SEM/岁':tt.SEM_batch_adjusted.map(lambda z:f'{z:.3f}'),'调整后MDC95/岁':tt.MDC95_batch_adjusted.map(lambda z:f'{z:.3f}')})
 cc=c.copy();cctab=pd.DataFrame({'类型':cc.celltype,'间隔/天':cc.stratum,'50次ICC中位数 [2.5–97.5百分位]':[f'{m:.3f} [{lo:.3f}, {hi:.3f}]' for m,lo,hi in zip(cc.ICC_median,cc.ICC_repeat_p025,cc.ICC_repeat_p975)]})
 aa=a[a.dataset.eq('AIDA_external')];atab=pd.DataFrame({'类型':aa.celltype,'变体':aa.variant,'N':aa.n_donors,'MAE/岁':aa.MAE.map(lambda z:f'{z:.2f}'),'Pearson r':aa.Pearson_r.map(lambda z:f'{z:.3f}'),'平均误差/岁':aa.bias_mean_error.map(lambda z:f'{z:+.2f}'),'偏移²/MSE':aa.bias_squared_fraction_MSE.map(lambda z:f'{100*z:.1f}%'),'主要成分':aa.dominant_AIDA_error_component.map({'fixed_offset':'固定偏移','scatter':'散布'})})
 # All three cohorts' requested accuracy is visible without a45-row table.
 precision=[]
 for (ct,var),g in a.groupby(['celltype','variant']):
  row={'类型':ct,'变体':var}
  for ds,name in [('OneK1K_test20','OneK留出'),('Allen_K0_full','Allen全68人'),('AIDA_external','AIDA')]:
   z=g[g.dataset.eq(ds)].iloc[0];row[name]=f'{z.MAE:.2f} / {z.Pearson_r:.3f} / {z.bias_mean_error:+.2f}'
  precision.append(row)
 text=f'''# 单细胞免疫时钟可靠性审计：第4轮

**F2结论：{dec['F2_result']}（{dec['status']}，内部完整队列扩展）。** e为本轮主要变体，c为次要变体；在读取本轮结果前，于2026-10-27 02:28:44 UTC（本地2026-10-26）写入活动台账。MONO和NK在90天与一年共4项目标比较的ΔICC95%CI均高于0，B的两个间隔均未确认。固定种子20261025。

{md(status[['step','status','conclusion']])}

## 恢复、模型范围与门槛

round3_checkpoint.zip SHA256为`01622a44f2a6e6f3c38a645f6ee90a1e87989de7b0d2efa9e37234632be2240a`，146个清单文件均核验。{link('F2_prespecified.json')}保留先写入的门槛与主要变体；{link('round4_execution_contract.json')}固定执行细节。

**T分析使用发表的sc-ImmuAging五类评分及RAM；F2使用第3轮OneK-only复刻a/e及新的交叉拟合c。它们不是同一个基线模型，ICC不能相互替换。** e与a的OneK权重完全冻结，68人只作评估，没有重新按S/E选择或拟合模型。旧S/E仅用于复用缓存：原E的165访视与补齐的155访视合为320访视、68人、4,609,769个五类细胞；native QC未剔除新补齐细胞。原E的1650个a/e预测复算最大差5.12×10⁻¹³年，支持输入与冻结模型一致。

a是OneK-only Pearson筛选＋LASSO分支的匹配复刻，不是作者完整五队列模型。a/c输入为log1p CP10K、100个15细胞pseudocell的连续均值；e先合并donor×visit×类型原始计数，再取log2(1+CPM)。**e同时改变聚合方式与表达尺度，F2不单独识别聚合操作本身的效果。** 所有分母均为完整释放raw library，未仅用选中基因重算。

F2要求B/MONO/NK至少一类在90天或一年，配对ΔICC的95%CI不含0；本轮全部满足项都指向改善。F2没有继承第3轮F的点增幅≥0.15要求，也没有新增精度门槛。本轮四个阳性目标CI为逐项95%CI，未按六项比较校正总体错误率。全68人包含此前E且主要变体是在第3轮之后预设，**这是内部扩展确认，不是独立外部复现**。

## T：技术与生物可靠性

第3轮31/42项批次占优全部是RAM：Down原始8项、年龄调整8项；Up年龄调整8项、原始7项。8种RAM亚群为CM CD4、CM CD8、Core naive CD4、Core naive CD8、GZMB−CD27+ EM CD4、GZMB−CD27− EM CD4、GZMK+CD27+ EM CD8及Naive CD4 Treg；Up原始评分中唯一不以批次为最大项的是GZMB−CD27+ EM CD4。**这42个方差模型合并全部320访视拟合，不带独立间隔标签。** {link('T_31_batch_dominant_inventory.csv')}逐行列明31项及其合并范围。

按冻结批次台账配对，同/跨批次为7–9天18/15人、75–105天5/23人、330–400天0/62人。全年没有同批次参照，90天同批次仅5人，记SMALL_N。比如原始钟B的90天残差ICC，同批次0.926、跨批次−0.099；这些是不同供者子集，不能将差值视为随机化批次效应。{link('T_pair_batch_ledger.csv')}和{link('Figure2_same_cross_batch_ICC.csv')}包含21个时钟/亚群系列的全部配对数与ICC，n<3保留NA而不填0。

采用用户允许的混合模型路线，在**时钟评分层面**拟合预测年龄~年龄+截距+donor随机效应+批次随机效应+访视残差，将A1条件细胞抽样/pseudocell方差作为已知对角噪声。REML估计方差后，逐donor用其余donor估计该批次BLUP，校正评分为原观测预测减去批次BLUP，保留访视残差，再重算LODO年龄残差、ICC、SEM与MDC。没有把评分替换成donor拟合均值，没有声称做了RNA计数矩阵ComBat-seq。304/320访视可由其他donor提供批次支持；其余16访视的未知批次效应按0处理。{link('T_batch_adjusted_predictions.csv')}同时保存原值、其他donor批次效应及仅补充用的同样本BLUP。

{md(techtab)}

所有三个间隔、五类的调整后指标及95%CI见{link('Figure3_batch_adjusted_reliability.csv')}。2000次donor bootstrap重新拟合LODO，**CI条件于已估计的批次校正，不传播批次模型估计误差**。以上方差份额来自本轮含固定年龄的新模型，不能反向替换第3轮原表。CD8短间隔ICC由0.922变为0.914，调整并非总会提高每个结果。

**技术贡献的可支持表述：** 模型将B/CD4/CD8/MONO/NK非donor方差中的32.7%/64.8%/46.8%/58.1%/58.7%归于批次关联部分。批次占总方差比例更小，不能混用分母。批次关联可能包含季节、时间与处理因素；访视残差混合自然波动和未建模技术误差。因此严格“纯技术/纯生物”拆分仍不可识别，不能把剩余百分比全部称为生物波动。B的90天调整ICC仍仅0.238，表明批次校正没有消除全部不稳定性。{link('T_batch_contribution_by_type_interval.csv')}另列模型预期批次配对差异方差/观测差异MSE，作为模型依赖的尺度描述。

![图2]({O/'Figure2_batch_structure.png'})
![图3]({O/'Figure3_batch_adjustment.png'})

## F2主要变体e：确认

使用全部68人，90天28对、一年62对。每次bootstrap从68个donor聚类重抽，a/e共享同一组multiplicity；每次重新拟合排除被评估donor全部访视的LODO年龄模型。直接取2000个ICC(e)−ICC(a)的2.5、97.5百分位，**没有用两个独立CI的端点相减**。全部20项ICC及10项配对差值均2000次有效。

{md(tab)}

非F2目标类型也完整报告：CD4两个间隔的差值CI均为正，CD8两个间隔均含0；它们不参与F2。{link('Figure4_F2_full68_a_e_ICC.csv')}给a/e自身ICC及CI，{link('Figure4_F2_paired_delta_ICC.csv')}给所有五类差值及配对CI。门槛判定见{link('F2_decision.json')}。

## 次要变体c：50次×2折交叉拟合

按年龄带随机分层，每次训练/评估各34人，保存全部3400个donor折ID。每个训练半队列只有12–21名7–9天配对donor；仅由这些训练donor计算候选基因pseudobulk log2(1+CPM) ICC。候选5331基因为OneK各训练折/完整训练年龄筛选排名的事前并集，不使用Allen表现选择候选。每个Allen半队列mask与OneK的折内Pearson排名取交集，固定第3轮a的基因数规则，在原OneK训练80%/五折中重调λ，再拟合。共500个模型，未使用对应Allen评估donor、OneK留出20%或AIDA结果调参。

每次重复，每个Allen donor的全部访视由相反半队列mask训练的模型预测；随后在全68人重新做LODO并计算ICC。以下为50次拆分的中位数与百分位，**不是donor bootstrap置信区间**。

{md(cctab)}

另保存每donor×visit平均50个折外预测后的ensemble ICC。外部精度中的c则为100个折模型的均值；它不是第3轮单个S筛选c。{link('c_50x2_donor_fold_ledger.csv')}、{link('c_crossfit_mask_QA.csv')}、{link('c_500_model_training_QA.csv')}、{link('c_Allen_50_repeat_OOF_predictions.csv')}和{link('c_OOF_ensemble_ICC.csv')}可逐层核查。五类均成功恢复旧c固定预测，最大差均小于0.001年。**c全程是次要分析，不用于决定F2。**

![图4]({O/'Figure4_repair_confirmation.png'})

## 三队列精度与AIDA偏倚

下表每个单元依次为 **MAE/岁、Pearson r、平均误差/岁**。Allen对donor等权，再平均其访视绝对误差；r以donor均值计算。OneK留出各类197人，MONO195人；Allen均68人；AIDA B619、CD4/CD8/MONO各589、NK618人。

{md(pd.DataFrame(precision))}

{link('Figure5_accuracy_bias_decomposition.csv')}保留全部45项数字及分母，{link('Figure5_accuracy_wide.csv')}提供宽表。AIDA 4044个潜在候选注释覆盖4043个，但唯一未注释C8orf44在a/c/e任何实际模型中权重均为0；**全部实际非零特征覆盖100%**。{link('AIDA_a_c_e_feature_coverage.csv')}逐模型族核验，不能将外部误差归因于这一个未用特征。

{md(atab)}

按恒等式MSE=平均误差²+中心化误差方差拆分。对五类等权平均MSE后，a/c/e的类型特异固定偏移份额分别为86.6%/62.4%/47.7%。a多数大误差由固定偏移主导；c混合；e的B/CD4/CD8以散布主导，MONO/NK偏移仍占54.0%/73.1%。这不是将五类不同方向偏倚先相加抵消。全局均值误差移除后的MAE仅作使用外部结局的诊断，**没有用于模型校准或F2判定**。平移不改变Pearson r，因此a的B、c的MONO/NK等低相关不能仅靠修正均值解决。

![图5]({O/'Figure5_accuracy_and_bias.png'})

## B：采样间隔核对

[sc-ImmuAging原文](https://www.nature.com/articles/s43587-025-00819-z)Methods的“Differentially expressed genes identification”明确T0至T3m；“BCG cohort”指向EGAS00001006990。[原始300BCG单细胞研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC10242447/)STAR Methods /300BCG study design同样描述接种前采血及三个月后追加采血。短原文摘录、章节及链接保存在{link('B_sampling_interval_primary_sources.csv')}，归档sc-ImmuAging正文的短语已用脚本复核。原始研究描述39人，本次源数据实际38对，未补造第39人。

因此匹配75–105天层：CD8年龄加速SEM **3.376807年 [2.381923,4.265491]**；第3轮Holm临界SEM **1.8年**，SEM/临界值 **1.876004**、CI下限/临界值 **1.323291**；反向临界值/SEM为0.533048。{link('B_SEM_interval_match.csv')}。7–9天只作敏感性，不能代替主时间尺度。公开Methods给的是约三个月，未提供38人的逐人准确采样日期，不能断言每人恰好90天；RPMI培养与HIHA处理也不同，时间尺度匹配并不证明外源SEM完全可迁移。本轮没有重跑BCG零假设模拟，也没有将未拒绝解释为BCG平均效应为零。

## 交付

图2–5各有PNG、SVG和作图CSV；{link('REPRODUCE.md')}给出脚本顺序和缓存边界，{link('manuscript_methods_results_zh.md')}更新方法与结果段落。versions_and_seed.json、input_files_SHA256.csv、SHA256SUMS.csv与checkpoint receipt绑定版本、输入和交付。大型原始计数不重复打包；500个次要c模型、实际使用的a/e模型、折外预测、配对bootstrap及基因ICC masks均随归档提供。旧round3/A1的结果、权重和归档保持不变。
''';(O/'analysis_report_zh.md').write_text(text)
 method=f'''# 方法与结果段落：第4轮更新

方法：在恢复并核验第3轮归档后，我们在任何本轮计算前指定donor pseudobulk变体e为主要修复方案，ICC≥0.7基因筛选变体c为次要方案。F2预设要求在完整68人K0队列的90天或一年间隔，B、MONO或NK至少一种细胞的ICC(e)−ICC(a)配对donor bootstrap95%CI不含0。冻结OneK训练的a/e权重不变，保留原28/62对稳态访视。采用2000次共享donor聚类重抽，对每次重抽分别重新拟合LODO年龄残差模型，直接从配对ICC差值分布取百分位CI。F2为六项逐项比较的内部扩展确认，未作总体多重检验保证，也不是新外部确认队列。

复刻范围：a为OneK-only Pearson年龄筛选＋LASSO分支，不等同于原作者五队列完整模型。a/c采用log1p CP10K、每donor×visit×类型100个15细胞pseudocell的连续均值；e采用合并原始计数的log2(1+CPM)。e同时改变聚合与表达尺度，不能将其差异仅归于去掉pseudocell。完整raw library用于所有归一化分母。

批次方法：原发表sc-ImmuAging评分采用含年龄固定效应、donor/批次交叉随机效应和访视残差的模型，条件细胞/pseudocell重采样方差作为已知对角噪声。利用其他donor估计每个待评估donor的批次BLUP，从原评分中减去该项后重算年龄预测和LODO可靠性，保留访视残差。无其他donor批次支持的16/320访视不作未知批次校正；批次校正CI条件于估计的校正参数。本方法是评分层面的随机效应敏感性分析，不是RNA矩阵ComBat-seq。非随机批次/时间结构及缺乏同批次一年配对，使严格技术/生物分离不可识别。

次要交叉拟合：进行50次年龄分层两折，每个Allen训练半队列的12–21名短间隔配对donor独立计算候选基因ICC。仅据该mask和固定的OneK训练数据建立c模型，共500个模型，Allen评估半队列不参与基因筛选。每次重复拼接全68人的折外预测后重算LODO ICC，报告50次中位数及2.5–97.5百分位，区别于bootstrapCI。另对每人每访视50个折外预测取均值形成Allen ensemble，外部c精度使用100个折模型的均值。

结果：原31/42项批次最大成分全部属于RAM，基于合并320访视，非间隔模型。7–9天、90天、一年同/跨批次对数分别为18/15、5/23、0/62。评分混合模型将B、CD4、CD8、MONO、NK非donor方差的32.7%、64.8%、46.8%、58.1%、58.7%归于批次关联。90天原发表钟ICC调整前后分别为B0.107→0.238、CD4 0.757→0.833、CD8 0.817→0.860、MONO0.406→0.557、NK0.642→0.820。剩余波动不可全部视为生物学。

主要修复结果：MONO90天ICC由a的0.179提高到e的0.644，Δ=0.465（95%CI0.140–0.739）；一年由0.176提高至0.638，Δ=0.461（0.163–0.718）。NK90天由0.698提高至0.935，Δ=0.237（0.099–0.551）；一年由0.755提高至0.879，Δ=0.124（0.047–0.300）。四项CI均为正，达到预设F2。B的两个间隔均未确认，且c不参与该判定。50次交叉拟合中，c的MONO90天ICC中位数0.657（拆分百分位0.266–0.924）、一年0.698（0.378–0.918）；NK分别0.918（0.858–0.948）和0.888（0.839–0.931）。这些范围反映拆分不稳定性，不能当作供者重抽CI。

精度与解释：按细胞类型等权MSE分解，AIDA的a/c/e误差中固定偏移平方份额为86.6%/62.4%/47.7%。e的MONO和NK仍有14.72/21.95年的MAE及+12.68/+21.16年的平均误差，Pearson r仅0.288/0.332。均值偏移校正不会改变相关，可靠性改善不代表准确年龄标定或通用临床效度。

采样间隔：sc-ImmuAging与原始300BCG Methods均描述接种前T0至三个月后T3m的scRNA采样，故75–105天SEM是主时间尺度匹配。CD8残差SEM3.377年，其95%CI下限2.382年，分别为第3轮1.8年临界值的1.876和1.323倍。逐人准确日历间隔未公开，BCG培养及HIHA流程也不同，因此仍需保留外源SEM可迁移的假设。未拒绝指定波动零假设不证明平均BCG效应为零。
''';(O/'manuscript_methods_results_zh.md').write_text(method)
 (O/'figure_captions_zh.md').write_text('''# 图注

图2：左为第3轮合并访视的年龄调整评分方差成分，不能解读为每个间隔的方差模型；右为同/跨批次的LODO残差ICC，标签给配对donor数。90天同批次仅5对，一年没有同批次配对，不显示替代数值。

图3：原发表sc-ImmuAging在五类×三个间隔的批次调整前后ICC、SEM及MDC95。评分保留访视残差，减去由其他donor估计的批次BLUP。调整后CI条件于已拟合校正，未传播其估计误差；批次关联不等于技术因果作用。

图4：左为主要e相对冻结复刻a的完整68人配对ΔICC及2000次donor bootstrap逐项95%CI，蓝色为F2预设目标类型。右为次要c50次两折OOF的ICC中位数及拆分百分位，灰点为a；右侧区间不是bootstrapCI。90天28对、一年62对。所有类型均展示，不只展示阳性比较。

图5：a/c/e在固定OneK留出、Allen全队列及AIDA的MAE、donor均值Pearson相关；底部按MSE恒等式分解AIDA固定偏移平方与中心化误差方差的份额。c为重复交叉拟合ensemble，a/e为冻结单模型。偏移移除仅为诊断，未进行使用目标结局的模型校准。
''')
 (O/'REPRODUCE.md').write_text('''# Round4 reproduction

Run from the project root with work/analysis_env/bin/python, seed20261025. Restore immutable round3 and round2_A1 outputs and their source caches; hashes/versions are bound in the manifests. Do not overwrite prior scientific outputs.

1. Record F2 and primary e/secondary c in gate_ledger.json before calculation. F2_prespecified.json retains this exact record; round4_execution_contract.json freezes implementation details.
2. Check round3 archive and all146 manifest files. T_31_batch_dominant_inventory.csv and T_pair_batch_ledger.csv are direct source-table transformations; original42 variance rows are pooled, not interval-specific.
3. analyze_technical_reliability.py computes stratified ICC and the allowed score-level random-effect adjustment. It does not run ComBat-seq. Original native clock and retrained OneK a/e are distinct analyses.
4. prepare_round4_inputs.py S prepares remaining old S-cache donors; extra fetches only missing candidate genes for the19 old E short-pair donors. Native full raw-library QC and stable source IDs are retained. For speed, prefetch_K0_S.py0 and1 may download disjoint even/odd chunks with per-chunk locks, then run prepare_round4_inputs.py S to aggregate cached counts. Superseded serial/interrupted download logs are retained. Full-row read followed by local feature projection avoids slow noncontiguous TileDB gene-range reads without changing selected features.
5. prepare_c_crossfit.py freezes50 age-stratified donor split ledgers and caches OneK-only sufficient statistics. build_c_masks.py computes geneICC using only opposite heldout-fold donors. fit_c_crossfit.py trains500 secondary models; optional R4_C_TYPES and R4_BLAS_THREADS enable disjoint cell-type processes. Recovery QA reproduces old c fixed predictions before accepting the new fits. No OneK heldout age, AIDA or Allen heldout outcome chooses genes/penalties.
6. evaluate_F2.py assembles full68 inputs and evaluates fixed a/e using2000 shared bootstrap draws, SeedSequence([20261025,4202]); paired differences are taken within each draw before percentiles. c models do not enter F2.
7. evaluate_c_and_accuracy.py builds true OOF predictions for each Allen repeat,50-repeat ICC summaries and ensembles, then computes all45 a/c/e accuracy/bias rows. audit_AIDA_coverage.py verifies actual nonzero-feature coverage. Repeated-split percentiles are not bootstrap CIs.
8. audit_BCG_interval.py verifies the archived sc-ImmuAging Methods phrase and binds the curated original300BCG primary Methods review. Exact individual collection dates are not available; nominal interval is3months. This step copies the round3 R1 SEM/critical ratios without rerunning the null simulation.
9. plot_round4.py, build_round4_report.py and package_round4.py produce figures, manuscript text, source preservation checks, versions, file hashes and checkpoint.

Large raw count and OneK statistic caches are not bundled. All500 c models, the10 fixed a/e models, ICC masks, OOF predictions and paired bootstrap distributions are bundled. Model prediction is intercept+sum(input*coefficient); use the expression_unit stored in each NPZ. PCA is not used in this round. c external accuracy is the100-model average; Allen c is50 true OOF predictions per donor/visit. All residual age calibrations exclude the held donor's complete visit history.
''');print('ROUND4_REPORT_COMPLETE',dec['F2_result'])
if __name__=='__main__':main()
