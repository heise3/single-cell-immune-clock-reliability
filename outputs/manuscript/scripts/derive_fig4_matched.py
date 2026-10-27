"""Derive within-person change for Allen 75-105 d pairs (no recorded vaccination between draws) and the
published BCG CD8 changes, used in Fig. 4a. Post hoc descriptive analysis added at manuscript stage."""
import os, pandas as pd, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); R=os.path.abspath(os.path.join(HERE,'..','..')); DER=os.path.join(R,'manuscript','derived')
x=pd.read_csv(f'{R}/round2_A1/main_predictions_with_LODO.csv'); x=x[x.clock=='sc-ImmuAging']
pp=x.set_index(['gsm','celltype']).prediction
led=pd.read_csv(f'{R}/round4/T_pair_batch_ledger.csv')
vac=pd.read_csv(f'{R}/manuscript_review_20261027/Allen_pair_vaccine_history_check.csv')
rows=[]
for st in ['75-105d','7-9d']:
    l=led[led.stratum==st].merge(vac[vac.stratum==st][['donor','any_recorded_vaccine_after_start_draw_before_end_draw']],on='donor',how='left')
    assert l.any_recorded_vaccine_after_start_draw_before_end_draw.notna().all()
    l=l[~l.any_recorded_vaccine_after_start_draw_before_end_draw.astype(bool)]
    for ct in ['B','CD4T','CD8T','MONO','NK']:
        for _,r in l.iterrows():
            b=pp[(r.gsm1,ct)]; f=pp[(r.gsm2,ct)]
            rows.append(dict(source='Allen',stratum=st,celltype=ct,donor=r.donor,interval_days=r.interval_days,batch_relation=r.batch_relation,baseline=b,followup=f,change=f-b))
bcg=pd.read_csv(f'{R}/round2_A1/BCG_38_paired_CD8.csv')
for _,r in bcg.iterrows():
    rows.append(dict(source='BCG',stratum='3 months post-BCG',celltype='CD8T',donor=r.donor_id,interval_days=np.nan,batch_relation='na',baseline=r.baseline,followup=r.followup,change=r.change))
D=pd.DataFrame(rows); D.to_csv(f'{DER}/Fig4a_matched_interval_change.csv',index=False)
S=D.groupby(['source','stratum','celltype']).apply(lambda g: pd.Series(dict(n=len(g),n_decline=int((g.change<0).sum()),pct_decline=100*(g.change<0).mean(),
    median_decline_among_decliners=g.change[g.change<0].median(),n_decline_ge_2p75=int((g.change<=-2.75).sum()),mean_change=g.change.mean(),median_change=g.change.median(),sd_change=g.change.std(ddof=1),
    r_baseline_change=np.corrcoef(g.baseline,g.change)[0,1]))).reset_index()
S.to_csv(f'{DER}/Fig4a_matched_interval_summary.csv',index=False)
print(S.round(2).to_string())
