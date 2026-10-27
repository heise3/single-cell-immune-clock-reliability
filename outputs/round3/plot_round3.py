#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[2];O=P/'outputs/round3'
def save(fig,name):
 fig.savefig(O/(name+'.svg'));fig.savefig(O/(name+'.png'),dpi=240);plt.close(fig)
def main():
 x=pd.read_csv(O/'R1_SEM_sensitivity_curve.csv');rr=pd.read_csv(O/'R1_critical_SEM_vs_CD8_estimates.csv');crit=json.loads((O/'R1_summary.json').read_text())['critical_SEM_Holm_grid'];fig,axs=plt.subplots(1,2,figsize=(12,4.7))
 for key,label,col in [('negative_count_Holm_p','Negative count','#0072B2'),('negative_median_Holm_p','Negative-subgroup median','#D55E00')]:axs[0].plot(x.SEM,x[key],label=label,color=col)
 axs[0].set_yscale('log');axs[0].axhline(.05,color='black',ls='--',lw=1);axs[0].axvline(crit,color='#777777',ls=':',label=f'Critical grid SEM = {crit:.1f}');axs[0].set_xlabel('Assumed SEM (years)');axs[0].set_ylabel('Two-sided empirical p after Holm(2)');axs[0].set_ylim(1e-5,1.2);axs[0].legend(frameon=False,fontsize=9);axs[0].set_title('BCG paired CD8: 100,000 simulations per SEM')
 rr['interval_order']=rr.stratum.map({'7-9d':0,'75-105d':1,'330-400d':2});rr['metric_order']=rr.metric.map({'prediction':0,'age_acceleration':1});rr=rr.sort_values(['metric_order','interval_order']);labels=[]
 for i,row in enumerate(rr.itertuples()):
  col='#0072B2' if row.metric=='age_acceleration' else '#999999';axs[1].hlines(i,row.SEM_CI_low,row.SEM_CI_high,color=col);axs[1].plot(row.SEM,i,'o',color=col);labels.append(('Age acceleration' if row.metric=='age_acceleration' else 'Raw prediction')+' / '+row.stratum)
 axs[1].set_yticks(np.arange(len(rr)),labels,fontsize=8);axs[1].invert_yaxis();axs[1].axvline(crit,color='black',ls='--',lw=1);axs[1].set_xlabel('Estimated CD8 SEM and donor-bootstrap 95% CI (years)');axs[1].set_title('All A1 CD8 SEM estimates versus critical grid value')
 for ax in axs:ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12)
 fig.tight_layout();save(fig,'Figure_R1_SEM_tipping_point')
 if (O/'R2_E_reliability_results.csv').exists():
  d=pd.read_csv(O/'R2_E_reliability_results.csv');types=['B','CD4T','CD8T','MONO','NK'];variants=['a','b','c','d','e'];colors=['#555555','#0072B2','#009E73','#CC79A7','#D55E00'];fig,axs=plt.subplots(2,5,figsize=(15,6),sharex=True,sharey=True)
  for row,st in enumerate(['75-105d','330-400d']):
   for col,ct in enumerate(types):
    ax=axs[row,col];z=d[d.celltype.eq(ct)&d.stratum.eq(st)].set_index('variant').reindex(variants)
    for i,var in enumerate(variants):
     r=z.loc[var];ax.hlines(i,r.ICC_CI_low,r.ICC_CI_high,color=colors[i]);ax.plot(r.ICC2_1,i,'o',color=colors[i],ms=4)
    if row==0 and ct in ['B','MONO','NK'] and np.isfinite(z.loc['a','ICC2_1']):ax.axvline(z.loc['a','ICC2_1']+.15,color='black',ls=':',lw=.8)
    ax.axvline(0,color='#BBBBBB',lw=.6);ax.set_title(f'{ct} / {st}\npaired donors n={int(z.iloc[0].n_paired_donors)}',fontsize=9);ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.12);ax.set_yticks(np.arange(5),variants);ax.set_xlim(-.6,1.1)
  axs[0,0].invert_yaxis();fig.supxlabel('LODO age-acceleration ICC(2,1), donor-bootstrap 95% CI');fig.supylabel('Variant (a–e)',fontsize=10);fig.suptitle('a baseline; b gene ICC≥0.5; c gene ICC≥0.7; d PCA; e donor pseudobulk\nDotted line: baseline +0.15 in prespecified F target types',fontsize=11);fig.tight_layout();save(fig,'Figure_R2_E_reliability')
  gate=pd.read_csv(O/'R2_F_gate_by_type.csv');gate=gate[gate.variant.ne('a')];fig,ax=plt.subplots(figsize=(7,5));palette={'B':'#0072B2','MONO':'#D55E00','NK':'#009E73'};markers={'b':'o','c':'s','d':'^','e':'D'}
  for r in gate.itertuples():
   ax.scatter(r.delta_ICC_vs_replicated_baseline,r.MAE_increase,color=palette[r.celltype],marker=markers[r.variant],s=50);ax.annotate(r.celltype+'/'+r.variant,(r.delta_ICC_vs_replicated_baseline,r.MAE_increase),xytext=(4,3),textcoords='offset points',fontsize=7)
  ax.axvline(.15,color='black',ls='--',lw=1);ax.axhline(1.5,color='black',ls='--',lw=1);ax.set_xlabel('90-day residual ICC increase versus matched baseline');ax.set_ylabel('OneK1K held-out MAE increase (years)');ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15);ax.set_title('Prespecified F: the same variant must qualify in at least two types');fig.tight_layout();save(fig,'Figure_R2_F_gate')
 print('PLOTS_COMPLETE')
if __name__=='__main__':main()
