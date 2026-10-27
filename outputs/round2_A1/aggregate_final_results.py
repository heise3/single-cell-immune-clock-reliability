#!/usr/bin/env python3
"""Run only after ordered steps2–5 have produced their required artifacts."""
from pathlib import Path
import json,datetime
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1'
def main():
 required=['reliability_results.csv','variance_decomposition.csv','main_predictions_with_LODO.csv','G1_prime_results.csv','RAM_runtime_gene_QA.csv','main_native_glmnet_validation.csv','PHARE_real_K0_run.json','siAge_alternate_routes.json','RTM_null_results.csv','RTM_summary.json','novelty_status.json','novelty_top5.csv','sensitivity87/reliability_results.csv']
 for name in required:assert (O/name).exists(),('Required ordered work remains',name)
 g1=pd.read_csv(O/'G1_prime_results.csv');r=pd.read_csv(O/'reliability_results.csv');v=pd.read_csv(O/'variance_decomposition.csv');pred=pd.read_csv(O/'main_predictions_with_LODO.csv');rtm=pd.read_csv(O/'RTM_null_results.csv');src=json.loads((O/'RTM_summary.json').read_text());nov=json.loads((O/'novelty_status.json').read_text());phare=json.loads((O/'PHARE_real_K0_run.json').read_text());si=json.loads((O/'siAge_alternate_routes.json').read_text())
 retained=set(g1.loc[g1.G1_prime_status.eq('PASS'),'celltype']);a=r[(r.clock=='sc-ImmuAging')&(r.metric=='age_acceleration')&r.celltype.isin(retained)];assert len(a)==len(retained)*3 and np.isfinite(a.ICC2_1).all();all_reliable=bool(a.ICC2_1.ge(.8).all());rr=rtm[rtm.role.eq('PRIMARY')];assert len(rr)==1;rr=rr.iloc[0];rejected=str(rr.null_rejected).lower()=='true';assert np.isfinite(rr.count_Holm_p) and np.isfinite(rr.negative_median_Holm_p)
 assert pred[pred.clock.eq('sc-ImmuAging')].donor.nunique()==68 and pred[pred.clock.eq('sc-ImmuAging')].gsm.nunique()==320
 assert r[r.clock.eq('RAM_Up')].celltype.nunique()==8 and phare['status']=='RUN_COMPLETED_REAL_MAIN_COHORT_SAMPLE'
 verdict='KILL' if all_reliable and rejected else 'PASS';legs=[]
 if not all_reliable:legs.append('可靠性低：至少一项保留类型/预定间隔的年龄加速ICC<0.8')
 if not rejected:legs.append('均值回归/波动零假设未被主检验拒绝')
 decision={'seed':20261025,'recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'K6_verdict':verdict,'all_retained_scImmuAging_age_acceleration_ICCs_at_least_0_8':all_reliable,'RTM_primary_null_rejected_after_Holm':rejected,'supporting_legs':legs,'primary_RTM_interval':'75-105d','retained_types':sorted(retained),'ICC_rows_evaluated':len(a),'rule':'unchanged user rule; point ICC threshold0.8 and any of two Holm-adjusted primary empirical p values<0.05','novelty_status':nov['status'],'novelty_complete_overlap':nov.get('complete_overlap'),'clock_completion':'sc-ImmuAging plus native-recipe RAM; PHARE one-sample execution, siAge routes separately reported'}
 if nov.get('complete_overlap') is True:decision['overall_status']='KILL_K2_FULL_NOVELTY_OVERLAP';decision['K6_interpretation']='K2 stop takes priority; retained statistical results are diagnostic only'
 elif nov['status']=='COMPLETE_NO_FULL_OVERLAP':decision['overall_status']=verdict
 else:decision['overall_status']='K6_'+verdict+'_NOVELTY_COVERAGE_LIMITED'
 (O/'decision.json').write_text(json.dumps(decision,indent=2,ensure_ascii=False))
 # Scope-aware comparison to the actual first-round pilot.
 old=pd.read_csv(P/'outputs/K4_reliability.csv').set_index('celltype');same=pd.read_csv(O/'same20_repair_and_rounding_comparison.csv');rows=[]
 for ct in sorted(retained):
  repaired=same[(same.celltype==ct)&same.analysis.eq('repaired_unrounded')].iloc[0];raw_old=same[(same.celltype==ct)&same.analysis.eq('old_unrounded')].iloc[0]
  for stratum in ['7-9d','75-105d','330-400d']:
   x=a[(a.celltype==ct)&a.stratum.eq(stratum)].iloc[0];q=r[(r.clock=='sc-ImmuAging')&(r.celltype==ct)&r.metric.eq('prediction')&r.stratum.eq(stratum)].iloc[0]
   rows.append({'celltype':ct,'stratum':stratum,'round1_donors':20,'round1_rounded_raw_ICC':old.loc[ct,'ICC_2_1'],'round1_rounded_MDC95':old.loc[ct,'MDC95'],'same20_old_unrounded_ICC':raw_old.ICC2_1,'same20_repaired_unrounded_ICC':repaired.ICC2_1,'same20_repaired_unrounded_MDC95':repaired.MDC95,'round1_missing_nonzero_genes':old.loc[ct,'missing_features_max'],'same20_repaired_missing_nonzero_genes':repaired.max_missing_nonzero_genes,'K0_primary_total_donors':68,'K0_interval_pair_donors':x.n_donors,'K0_raw_ICC':q.ICC2_1,'K0_age_acceleration_ICC':x.ICC2_1,'K0_age_acceleration_SEM':x.SEM,'K0_age_acceleration_MDC95':x.MDC95,'scope_note':'Only same20 columns isolate repair/rounding. K0 changes population, released labels, gene universe and estimand; cross-population differences are not attributed to gene repair alone.'})
 pd.DataFrame(rows).to_csv(O/'round1_to_A1_change_table.csv',index=False)
 eff=json.loads((P/'outputs/published_effect_magnitudes.json').read_text());effectrows=[]
 for ct,effects in [('CD8T',[('BCG negative-change subgroup median',abs(eff['BCG']['AR_median_years']))]),('MONO',[('COVID mild-control, exact age63 subgroup',eff['COVID_MONO']['mild_minus_control_years']),('COVID severe-control, exact age63 subgroup',eff['COVID_MONO']['severe_minus_control_years'])])]:
  for x in a[a.celltype.eq(ct)].itertuples():
   for label,val in effects:effectrows.append({'celltype':ct,'stratum':x.stratum,'MDC95_individual':x.MDC95,'MDC95_CI_low':x.MDC95_CI_low,'MDC95_CI_high':x.MDC95_CI_high,'published_source_statistic':label,'effect_magnitude_years':val,'MDC_over_effect':x.MDC95/val,'interpretation':'Descriptive scale comparison only: individual MDC and selected subgroup/group contrast are different estimands; not a test of a population treatment effect.'})
 pd.DataFrame(effectrows).to_csv(O/'MDC_vs_published_effect_scale.csv',index=False)
 # Shared raw-age predictions must be identical despite separate age calibrations.
 sens=pd.read_csv(O/'sensitivity87/main_predictions_with_LODO.csv');share=pred[pred.clock.eq('sc-ImmuAging')].merge(sens[sens.clock.eq('sc-ImmuAging')],on=['clock','celltype','gsm'],suffixes=('_K0','_87'),validate='one_to_one');share['raw_prediction_abs_difference']=abs(share.prediction_K0-share.prediction_87);share[['celltype','gsm','raw_prediction_abs_difference']].to_csv(O/'shared_visit_prediction_invariance_QA.csv',index=False);assert share.raw_prediction_abs_difference.max()<1e-8
 statuses=[{'step':'A1','status':'PASS','evidence':'A1_amendment.json; gate_ledger.json','conclusion':'修正先于计算，旧5人检验保留并取代'}, {'step':'0','status':'PASS','evidence':'68 donors,320 visits,90 batches; full observation_joinid and frozen metadata identity checks','conclusion':'主分析始终使用冻结K0；87人口径独立敏感性'}, {'step':"G1 prime",'status':'PASS' if g1.G1_prime_status.eq('PASS').all() else 'TYPE_SPECIFIC_EXCLUSION','evidence':'G1_prime_results.csv;200 donors/type,ages19–97','conclusion':'按r与方差比放行；偏倚只记录'}, {'step':'2','status':'PASS','evidence':'reliability_results.csv;variance_decomposition.csv;donor bootstrap2000','conclusion':'预测年龄与LODO年龄加速均完成三间隔分析'}, {'step':'3','status':'PASS','evidence':'sc-ImmuAging+RAM completed;PHARE real-sample runtime;siAge_alternate_routes.json','conclusion':'完成A1规定的两个可靠性完成项；其他时钟状态单独披露'}, {'step':'4','status':'PASS_EXECUTED','evidence':f'90d Holm p(count)={rr.count_Holm_p:.6g};p(negative-median)={rr.negative_median_Holm_p:.6g}','conclusion':'主检验拒绝指定零假设' if rejected else '指定零假设仍与数据相容'}, {'step':'5','status':nov['status'],'evidence':'novelty_top5.csv;novelty_search_log.csv','conclusion':'完全重叠' if nov.get('complete_overlap') is True else '成功读取的结果中未检出完全重叠；覆盖缺口按日志披露'}, {'step':'6','status':decision['overall_status'],'evidence':'decision.json','conclusion':'；'.join(legs) if legs else '两条KILL条件同时成立'}];pd.DataFrame(statuses).to_csv(O/'step_status.csv',index=False)
 print(json.dumps(decision,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
