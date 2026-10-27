#!/usr/bin/env python3
"""Render evidence tables and manuscript prose from frozen round5 result files."""
from pathlib import Path
import json,pandas as pd,numpy as np
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';S={'75-105d':'约90天','330-400d':'约1年'};CTS=['B','CD4T','CD8T','MONO','NK']
def md(d):
 d=d.copy().fillna('NA');return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v).replace('|','/') for v in r)+' |' for r in d.itertuples(index=False,name=None))
def read(n):return pd.read_csv(O/(n+'.csv'))
def main():
 p=read('P1_e_vs_published_paired_delta');h=read('P1_F2_Holm6');r2=read('P2_covariate_between_donor_R2');age=read('P2_chronological_age_signal');adj=read('P2_CMV_sex_adjusted_ICC');ret=read('P2_retained_reliability_gain');vc=read('P3_three_model_variance_components');cov=read('P4_model_feature_coverage');search=json.loads((O/'P4_search_status.json').read_text());select=json.loads((O/'P4_external_selection_contract.json').read_text());ext=read('P4_external_technical_ICC');erun=json.loads((O/'P4_external_run_status.json').read_text())
 table=[]
 for r in p.itertuples():table.append({'细胞':r.celltype,'间隔':S[r.stratum],'配对donor':28 if r.stratum=='75-105d' else 62,'已发表':f'{r.ICC_published:.3f}','复刻a':f'{r.ICC_a:.3f}','e':f'{r.ICC_e:.3f}','e−已发表 [95%CI]':f'{r.delta_e_minus_published:+.3f} [{r.delta_CI_low:+.3f}, {r.delta_CI_high:+.3f}]'})
 t1=md(pd.DataFrame(table));hh=h[['celltype','stratum','delta_e_minus_a','bootstrap_tail_p','Holm6_p','passes_Holm6','null_centered_p_sensitivity','centered_Holm6_p_sensitivity','passes_centered_Holm6_sensitivity']].copy();hh.stratum=hh.stratum.map(S);hh=hh.round(6)
 ar=age.pivot(index='celltype',columns='model',values='Pearson_r').loc[CTS,['published','a','e']].reset_index().round(3);rr=r2[['celltype','model','CMV_R2','Sex_R2','joint_R2','joint_adjusted_R2']].round(4);aa=adj[['celltype','model','stratum','n_pairs','ICC_before','ICC_after_CMV_sex','adjusted_ICC_CI_low','adjusted_ICC_CI_high']].copy();aa.stratum=aa.stratum.map(S);aa=aa.round(4)
 vv=vc[vc.metric.eq('age_adjusted')].copy();components=['donor','visit_nested_plus_unmodelled_residual','batch','cell_sampling','pseudocell_algorithm'];vabs=vv[['celltype','model']+[c+'_variance' for c in components]].round(4);vfrac=vv[['celltype','model']+[c+'_fraction' for c in components]].round(4)
 et=ext[['celltype','clock','metric','n_donors','ICC2_1','ICC_CI_low','ICC_CI_high','bootstrap_valid']].round(4)
 status=[{'步骤':'P1','状态':'PASS（计算完成）','证据':'68人；28/62对；2000次配对重抽，a/e的20组bootstrap数组与第4轮一致','一句话结论':'对已发表模型仅NK两个间隔的逐项差值CI全为正；e对a主Holm通过4项，中心化敏感性仅MONO2项。'}, {'步骤':'P2','状态':'PASS（检查完成；年龄效度受限）','证据':'CMV+/−28/40；女/男38/30；68人协变量完整且跨访视不变；BMI缺失','一句话结论':'MONO/NK增益调整后保留93%–110%，但MONO e的真实年龄r=0.169，不能把可靠性直接解释为年龄效度。'}, {'步骤':'P3','状态':'PASS（分解完成）','证据':'30项=3模型×5类型×原始/残差；100细胞重采样和100a算法实现；原发表10项原样复用','一句话结论':'五类最大的非donor方差压缩均为批次；细胞抽样方差反而增加。'}, {'步骤':'P4','状态':'PASS（技术复现）；生物纵向NOT_RUN','证据':f'四库筛查；AIDA431人符合双库五类QC；预抽20人、{select["n_cells"]:,}细胞；20个ICC及CI','一句话结论':'合格技术重复已运行；e五类ICC点估计均低于原发表，普遍优势未复现；同次采血不能外推生物纵向。'}];pd.DataFrame(status).to_csv(O/'step_status.csv',index=False)
 report=f'''# 单细胞免疫时钟可靠性审计：第5轮

从 round4_checkpoint 恢复，固定种子 **20261025**。本轮合同于 **2026-10-27T08:46:52.274855+00:00** 写入 gate_ledger.json，早于本轮结果计算；AIDA20人子集于 **{select['recorded_UTC']}** 登记，早于拆分文库的预测。未增加模型训练，未改K0人群、阈值、已有权重或历史F2决定。恢复时逐一验证1124个第4轮文件，归档SHA256为 `258655b4f4a7f81522c3567a25fbaa66f5d85b2020bc3a629b46efd322d5ff1c`。

{md(pd.DataFrame(status))}

本轮收紧的是可支持的论断：e对复刻a的改善不能整体写成优于原发表时钟；MONO的年龄信号仍弱；NK的多重校正推断有bootstrap算法敏感性；已完成的是同次采血技术重复，其中e五类ICC点估计均低于原发表模型，未复现普遍优势，生物纵向外部复现仍未运行。这些限制不倒写第4轮既定F2记录，也不新增事后PASS阈值。

**预注册标签另行声明：用户询问的P2验证性/探索性标签仍为“待声明”。** 第3轮已在14对短间隔S集上计算过31,750个CD4基因的ICC及MSR/MSC/MSE中间方差量；第4轮50×2折又计算并保存100×5,331个CD4基因ICC和mask。因此不能声明此前未接触CD4基因层面方差信息。目前未查到逐CD4基因五项REML的既有执行记录，但这不消除ICC中间量的既有接触。本轮P2“PASS”仅指完成指定计算，不等于预注册验证性。证据和范围另见 `CD4_gene_variance_prior_exposure_declaration.md`。

## P1：已发表 / 复刻a / e 与多重校正

主队列仍为68名donor、320次稳态访视。90天采用75–105天窗口28对，一年采用330–400天窗口62对。沿用第4轮同一组68人聚类重抽2000次的多项分布权重（SeedSequence[20261025,4202]）。每次均重新拟合LODO，排除被留出原始donor的全部访视和全部bootstrap副本，从同次重抽的ICC差值直接取百分位CI。未用两个独立CI相减。所有a/e原始bootstrap数组在1e−11内重现。

{t1}

以上差值CI为逐项95%CI；对e−已发表未另规定全家族显著性门槛。B和MONO两个间隔的差值CI均跨0，只有NK两个间隔为正。a是OneK-only Pearson年龄筛选＋LASSO复刻，原发表模型来自作者五队列模型，两者必须区别。

第4轮e−a的六项组成一个Holm家族。预设主要p值是两侧percentile-bootstrap尾概率：`min(1,2*min((N(delta_boot<=0)+1)/2001,(N(delta_boot>=0)+1)/2001))`。这是一种近似bootstrap尾概率，不是严格零假设随机化检验。另按预登记报告中心化零假设敏感性p：`(N(abs(delta_boot−delta_point)>=abs(delta_point))+1)/2001`，单独Holm6。

{md(hh)}

主要方法通过MONO两项和NK两项；中心化敏感性只有MONO两项通过。NK的差值分布不对称，两种构造p值差异明显。不能把4项通过表述为对p值构造稳健的4项总体确认。

## P2：稳定协变量与年龄效度

直接解析归档GEO [GSE271896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271896) family SOFT。68人CMV与性别均完整且在稳态访视内不变：CMV阴性40、阳性28；女性38、男性30。全部已释放characteristic键中没有BMI、体重、身高，不插补BMI。

每人LODO年龄残差跨访视取等权均值后估计描述性R²，三种模型各五类：

{md(rr)}

CMV与性别合计R²不能视作两个单独R²之和；同时提供调整R²、unique增量和partial R²。e的NK单独性别R²=32.2%，CMV R²=4.5%，合计38.9%（调整R²=37.0%）。这是稳定协变量关联，不能据此判定其生物学意义无关衰老。

调整用嵌套LODO：对每个目标donor，年龄校准与CMV+性别回归均完全排除目标donor；训练供者的残差及协变量效应在这个目标donor排除的世界里构造。2000次bootstrap同时重新拟合这两层；保留所有68人和同一组配对。

{md(aa)}

MONO90天/一年e相对a的增益分别保留95.0%/104.7%；NK分别保留110.4%/93.3%。比例按调整后差值/调整前正差值计算、不截断100%；原始差值非正的B90天、CD8两间隔不给有利的“保留率”。去除协变量未消除这些ICC增益，仍不能证明剩余稳定差异都是衰老。

全Allen68人各模型原始预测年龄与真实年龄的Pearson r（每人取预测和真实年龄均值，独立N=68）：

{md(ar)}

e保留明显CD4/CD8年龄信号并改善B的年龄相关性。MONO e r=0.169（p=0.169）仍弱，NK e r=0.325（p=0.00679）仅中等；不能作“所有细胞类型均保留充分年龄效度”的结论，也不能把高残差ICC等同于衰老构念效度。

## P3：修复机制的方差分解

沿用A1/第3轮的条件REML协方差模型：donor、批次交叉随机项、访视嵌套＋无法另分的技术残差；逐donor×visit×celltype的细胞重抽样与算法实现方差作为已知对角噪声。原始预测和LODO残差共30项。已发表模型10行直接复用A1；a/e新增条件细胞bootstrap100次，a固定细胞集下重复100次100×15 pseudocell实现，e为确定性pseudobulk，算法项结构性为0。新细胞抽样对子模型a/e采用同一重抽细胞；e每次重合并原始计数和完整library总量。

年龄加速的绝对方差（年²；原始预测另见CSV）：

{md(vabs)}

相应占比：

{md(vfrac)}

e相对a五类最大的非donor方差压缩均为批次。B、MONO、NK的批次方差从40.670→2.582、10.891→0.568、17.207→1.764年²，下降93.7%、94.8%、89.7%；其原发表对应值为3.135、12.033、10.278年²，故B相对原发表的批次改善仅约17.6%，不能用相对a的94%替代。CD8 e的批次方差7.644高于原发表4.923。e的条件细胞抽样方差在五类均高于a，所以改善不来自细胞抽样项压缩。e同时改变聚合方式和表达尺度，并重新用OneK训练，不是孤立改变pseudocell的机制实验；观察性批次结构也不允许将访视残差直接称为纯生物变异。

## P4：外部重复数据检索与一次技术复现

记录的官方检索包括GEO 89个series题录/摘要、CELLxGENE {search['CELLxGENE_collections_screened']}个collection元数据、HCA全部86个blood项目，以及ArrayExpress明确AND查询的29个healthy和8个repeated命中。题录筛查与逐候选深核验分别保存；不是每个命中都完成原始文件结构审计。初始HCA size100的400、宽泛BioStudies OR检索、EPMC503与10x直连429都保留在receipt；HCA改size50并读两页blood命中，BioStudies改AND，未把访问错误记成数据不存在。

逐候选编号、N、间隔、平台、签字/DUA、矩阵与元数据状态见 `P4_candidate_eligibility.csv`。最接近者包括PALMO/GSE190992（4名donor、每周重复，不足10）、GSE275071（4名稳态donor×4GEM，32样本不是32人）、GSE247277（11人只有一次baseline scRNA）、GSE206265（scRNA0/1/28天，第二次已接种）、GSE304704（3人冻存技术比较）和商业kit比较（1人）。GSE158769有12份技术重复但仅分选memory T，不能作为B/MONO/NK全PBMC复现。

HCA [Terekhova健康血液图谱](https://explore.data.humancellatlas.org/projects/896f377c-8e88-463e-82b0-b2a5409d6fe4)是另一条强纵向线索：166人、317样本、3年收集期，HCA列出分类型处理矩阵，allPBMC tar约17.1GB；[AWS登记](https://registry.opendata.aws/singlecellhumanbloodatlas/)对Synapse路线标有controlled access。该候选的逐人年度间隔、原始整数count层和确切成对N本轮未核验，不把317样本当317个独立donor，也不将它判为无数据。按本轮只运行一个小数据源的范围，选择已有计数缓存且元数据可核实的AIDA技术重复。

[AIDA作者资源](https://prabhakarlab.github.io/AIDA.html)与[官方10x示例](https://www.10xgenomics.com/support/jp/software/loupe-browser/latest/tutorials/lb-sc-analyzing-public-data)说明每样本具有两份文库；来源逐细胞sample/library/suspension与donor标识连接到已固定CELLxGENE v2和Census2025-11-08。这里是同次采血、同PBMC pool、独立GEM GEX文库，通常在同一NovaSeq lane测序。它检验技术重测，不覆盖访视、采血及跨批次的生物变动。

在结果计算前，从431名满足两份文库且五类各≥30个QC细胞的donor中，以SeedSequence[20261025,5401]均匀抽20名；年龄21–66岁，40个donor×library单元、36,805细胞。每人两库按UUID排序作为重复标签，不代表统一实验员或统一日历时间。保持既有native QC、完整raw library分母和固定模型；只补取1145个已测但未在旧4044基因缓存内的特征，绝不把工程缓存缺失当生物零。e非零权重覆盖100%；原发表B/CD4/CD8/MONO/NK映射覆盖分别98.83%/99.14%/99.84%/98.60%/99.53%，残余采用作者native补零并逐基因列出，因此外部结果特指该注释和缺失策略下的模型应用。

技术重复ICC（20对；2000次donor bootstrap，LODO每次重拟合）：

{md(et)}

五类细胞的原始年龄和年龄加速ICC点估计，e均低于原发表模型。因此这次技术重复没有复现e普遍更可靠的优势；本轮未预设外部差值门槛、也未计算外部配对差值CI，不能将点估计差异称为显著劣效。它与P3“e降低批次方差但增加细胞抽样方差”的权衡相一致，不过并非该机制的因果验证。

原始年龄、年龄加速均报告CI；`age_acceleration`为LODO残差。AIDA在前轮用于固定模型的外部精度报告，从未训练或按其结果选模型；本轮拆分library的可靠性为新增检验。它不是全新的年龄精度测试集，也不补足90天/一年生物纵向独立确认。本轮已找到并运行合格的公开技术重复数据，因此不写“无合格公开复现数据”。

## 图表、重现与状态

图P1为三模型ICC、e−已发表CI及e−a两种Holm p；图P2为协变量R²、年龄r和调整后增益；图P3为方差占比与绝对批次方差；图P4仅表示AIDA技术重复。作图CSV和PNG/SVG/PDF均在本目录，SVG保留可编辑文字。

脚本：`extract_GEO_covariates.py`、`analyze_P1_P2.py`、`analyze_P3.py`、`search_P4*.py`、`check_AIDA_repeats.py`、`prepare_P4_technical.py`、`analyze_P4_technical.py`、`summarize_P4_search.py`、`plot_round5.py`、`build_report.py`。`README_reproduce.md`列出输入、运行顺序、seed子流、版本与哈希。本轮未训练新模型；a/e冻存权重与原发表系数及第4轮归档分别SHA绑定，输入大矩阵仍按manifest定位，不在新checkpoint复制。
'''
 (O/'round5_report_zh.md').write_text(report)
 methods='''# 方法与结果：第5轮更新草稿

方法：在冻结68名供者、320次稳态访视以及既有模型权重的前提下，补充比较原发表sc-ImmuAging、OneK-only复刻a和pseudobulk变体e。75–105天和330–400天分别保留28对、62对。采用相同的2000次供者聚类bootstrap，在每次重抽中重新拟合留一供者年龄校准，剔除被留出原始供者的全部副本。以配对差值分布计算ICC(e)−ICC(原发表)的百分位95%CI。对预定六项e−a比较，按预先登记的两侧bootstrap尾概率做Holm校正，并独立报告中心化零假设bootstrap p值的Holm敏感性分析，避免以计算方便或有利结果选择p值定义。复刻a不等于作者五队列模型；e同时改变聚合和表达尺度。

结果：e对原发表模型的残差ICC增益，B在90天/一年为0.349（95%CI −0.067至0.727）/0.138（−0.259至0.480），MONO为0.238（−0.072至0.510）/0.081（−0.165至0.289），NK为0.293（0.095至0.692）/0.170（0.067至0.322）。只有NK两个逐项CI完全为正。相对复刻a，主要bootstrap尾概率Holm校正通过MONO两项和NK两项；中心化敏感性仅MONO两项通过，NK两项的校正p均约0.312，提示其显著性依赖bootstrap p值构造。因此第4轮的内部F2确认记录保留，但不概括为对所有比较方法稳健地优于原发表模型。

方法：从GEO原始SOFT解析CMV和性别，保持每个供者总权重一致。以各供者跨访视LODO残差均值估计CMV、性别单独及联合R²。调整可靠性时采用嵌套留一供者：年龄校准和协变量回归均排除目标供者；在每次bootstrap重新拟合，防止目标供者经训练残差间接泄漏。GEO未释放BMI、体重或身高，未进行插补。

结果：68人中CMV阳性28人、阴性40人，女性38人、男性30人。e的NK残差供者间方差由CMV与性别合计解释38.9%，其中性别单独解释32.2%。调整后MONO的90天/一年增益保留95.0%/104.7%，NK保留110.4%/93.3%，表明这两个稳定协变量不足以解释ICC增益。相应地，e的年龄相关r在B、CD4、CD8、MONO、NK分别为0.605、0.888、0.828、0.169、0.325。MONO年龄效度仍弱；稳定残差不能直接等同于衰老信号。

方法与结果：沿用前轮条件REML模型，分解donor、批次、访视加未建模技术残差、细胞抽样和pseudocell算法五项。每个donor×visit×细胞类型进行100次细胞重抽样；a再进行100次固定细胞下的100×15 pseudocell随机实现，e算法项按确定性构造为0。e相对a五类的最大非donor方差压缩均来自批次；B、MONO、NK的年龄残差批次方差下降93.7%、94.8%、89.7%，但其细胞抽样方差升高。相对原发表模型，B批次下降仅17.6%，CD8批次方差反而较高。该条件观察性分解不支持把剩余访视项称为纯生物变异，也不隔离聚合方式与表达尺度的因果贡献。

'''
 exrows=[]
 for ct in CTS:
  q=ext[ext.celltype.eq(ct)&ext.metric.eq('age_acceleration')].set_index('clock');exrows.append(f"{ct}：原发表{q.loc['published','ICC2_1']:.3f}（{q.loc['published','ICC_CI_low']:.3f}–{q.loc['published','ICC_CI_high']:.3f}），e {q.loc['e','ICC2_1']:.3f}（{q.loc['e','ICC_CI_low']:.3f}–{q.loc['e','ICC_CI_high']:.3f}）")
 methods+='外部技术复现：在GEO、CELLxGENE、ArrayExpress和HCA完成有记录的筛查后，选取AIDA固定版本中20名健康成人的两份独立GEM GEX文库，共36,805细胞。子集按元数据及每类每库≥30细胞门槛预先抽取，年龄21–66岁；不训练或选择新模型。以相同native预处理和2000次供者bootstrap计算原始年龄及LODO年龄加速ICC。年龄加速结果为'+ '；'.join(exrows)+'。e在五类的ICC点估计均低于原发表，未复现其普遍可靠性优势；本轮未作外部差值显著性推断。原发表模型映射覆盖98.6%–99.8%，e覆盖100%，残余缺失沿用native零填策略。AIDA技术重复来自同次采血和同PBMC池、通常同一测序lane，因此该结果仅支持技术重复层面的外部检验；90天或一年生物学重测仍未完成。HCA健康血液图谱166人/317样本是后续可核验的纵向候选，本轮未重建其时间和整数count层。\n\n稿件结论措辞：本研究支持细胞类型、时间间隔、模型来源和技术环节共同决定免疫年龄指标的重测表现。pseudobulk方案的内部改善主要伴随批次方差压缩，但年龄构念效度和外部生物纵向泛化尚不充分，故应将“可靠性改善”与“衰老测量得到修复”分开表述。\n'
 methods+='\n预注册接触声明：预注册P2验证性/探索性标签仍为待声明。此前已计算CD4基因ICC及MSR/MSC/MSE中间方差量（第3轮31,750基因/14对；第4轮100折×5,331基因/每折12–21对），不能作未见数据声明。未查到逐基因五项REML的既有执行记录；该区别不消除先前接触。详见独立声明文件。\n';(O/'manuscript_methods_results_zh.md').write_text(methods);(O/'round5_decision.json').write_text(json.dumps({'round5_status':'COMPLETE_TARGETED_ANALYSES','historical_F2':'unchanged_confirmed_under_original_pointwise_rule','P1_Holm_primary_passes':h[h.passes_Holm6][['celltype','stratum']].to_dict('records'),'P1_centered_sensitivity_passes':h[h.passes_centered_Holm6_sensitivity][['celltype','stratum']].to_dict('records'),'superiority_over_published_pointwise_CI':['NK75-105d','NK330-400d'],'MONO_age_construct_validity':'limited_r0.1686_p0.169','P4':'AIDAtechnical_external_run_complete;allfive e ICC point estimates below published;biological_longitudinal_NOT_RUN','preregistration_P2_label':'待声明','new_training':False,'no_new_posthoc_gate':True},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
