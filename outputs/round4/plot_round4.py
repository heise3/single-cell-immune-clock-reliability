#!/usr/bin/env python3
from pathlib import Path
import shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[2];O=P/'outputs/round4';TYPES=['B','CD4T','CD8T','MONO','NK'];STS=['7-9d','75-105d','330-400d'];plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
def save(fig,name):
 fig.savefig(O/(name+'.svg'),bbox_inches='tight');fig.savefig(O/(name+'.png'),dpi=240,bbox_inches='tight');plt.close(fig)
def clean(ax):ax.spines[['top','right']].set_visible(False)
def main():
 v=pd.read_csv(P/'outputs/round3/Table2_variance_components_42.csv');v.to_csv(O/'Figure2_variance_components.csv',index=False);z=v[v.metric.eq('age_adjusted')].reset_index(drop=True);fig,axs=plt.subplots(1,2,figsize=(17,9),gridspec_kw={'width_ratios':[1.2,1]});left=np.zeros(len(z));components=[('donor','#4D4D4D','Donor'),('batch','#D55E00','Batch'),('visit_nested_plus_unmodelled_residual','#999999','Visit + residual'),('cell_sampling','#0072B2','Cell sampling'),('pseudocell_algorithm','#009E73','Pseudocell')]
 for c,col,label in components:axs[0].barh(np.arange(len(z)),z[c+'_fraction'],left=left,color=col,label=label);left+=z[c+'_fraction']
 axs[0].set_yticks(np.arange(len(z)),[r.clock.replace('RAM_Down','RAM Down').replace('RAM_Up','RAM Up')+' / '+r.celltype for r in z.itertuples()],fontsize=8);axs[0].invert_yaxis();axs[0].set_xlim(0,1);axs[0].set_xlabel('Fraction of variance (age-adjusted score)');axs[0].set_title('Pooled 320 visits: components are not interval-specific');axs[0].legend(frameon=False,ncol=3,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.06));q=pd.read_csv(O/'Figure2_same_cross_batch_ICC.csv');q=q[q.clock.eq('sc-ImmuAging')];ys=[];labs=[]
 for ct in TYPES:
  for st in STS:
   y=len(labs);ss=q[q.celltype.eq(ct)&q.stratum.eq(st)].set_index('batch_relation');labs.append(f'{ct} / {st}  ({int(ss.loc["same_batch","n_pairs"])}/{int(ss.loc["cross_batch","n_pairs"])} pairs)')
   for rel,col,offset in [('same_batch','#0072B2',-.12),('cross_batch','#D55E00',.12)]:axs[1].plot(ss.loc[rel,'ICC2_1'],y+offset,'o',color=col,ms=5,label=rel if y==0 else None)
 axs[1].set_yticks(np.arange(len(labs)),labs,fontsize=8);axs[1].invert_yaxis();axs[1].axvline(0,color='#BBBBBB',lw=.7);axs[1].set_xlim(-.2,1.05);axs[1].set_xlabel('LODO age-acceleration ICC(2,1)');axs[1].set_title('Same / cross batch pairs; 90-day same-batch n = 5');axs[1].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.06),ncol=2);[clean(ax) for ax in axs];fig.tight_layout(w_pad=3);save(fig,'Figure2_batch_structure')
 old=pd.read_csv(P/'outputs/round2_A1/reliability_results.csv');old=old[old.clock.eq('sc-ImmuAging')&old.metric.eq('age_acceleration')];new=pd.read_csv(O/'Figure3_batch_adjusted_reliability.csv');new=new[new.metric.eq('age_acceleration')];data=old.merge(new,on=['celltype','stratum'],suffixes=['_before','_adjusted'],validate='one_to_one');data.to_csv(O/'Figure3_plot_data.csv',index=False);fig,axs=plt.subplots(1,3,figsize=(14,7),sharey=True);labs=[]
 for ct in TYPES:
  for st in STS:
   row=data[data.celltype.eq(ct)&data.stratum.eq(st)].iloc[0];i=len(labs);labs.append(ct+' / '+st)
   for ax,key,lo,hi in zip(axs,['ICC2_1','SEM','MDC95'],['ICC_CI_low','SEM_CI_low','MDC95_CI_low'],['ICC_CI_high','SEM_CI_high','MDC95_CI_high']):
    ax.plot([row[key+'_before'],row[key+'_adjusted']],[i,i],color='#BBBBBB',lw=1)
    for tag,col,off in [('before','#777777',-.1),('adjusted','#0072B2',.1)]:ax.hlines(i+off,row[lo+'_'+tag],row[hi+'_'+tag],color=col,lw=.8);ax.plot(row[key+'_'+tag],i+off,'o',color=col,ms=4,label=tag if i==0 else None)
 axs[0].set_yticks(np.arange(len(labs)),labs);axs[0].invert_yaxis()
 for ax,label in zip(axs,['ICC(2,1)','SEM (years)','MDC95 (years)']):ax.set_xlabel(label);ax.grid(axis='x',alpha=.13);clean(ax)
 axs[0].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.1),ncol=2);fig.suptitle('Native sc-ImmuAging: subtract other-donor batch BLUP, retaining visit residual\nAdjusted CIs condition on fitted batch correction; technical causation is not identified',fontsize=11);fig.tight_layout();save(fig,'Figure3_batch_adjustment')
 if (O/'Figure4_F2_paired_delta_ICC.csv').exists() and (O/'Figure4_c_repeated_crossfit_summary.csv').exists():
  d=pd.read_csv(O/'Figure4_F2_paired_delta_ICC.csv');c=pd.read_csv(O/'Figure4_c_repeated_crossfit_summary.csv');d=d.merge(c,on=['celltype','stratum'],validate='one_to_one');d.to_csv(O/'Figure4_plot_data.csv',index=False);fig,axs=plt.subplots(1,2,figsize=(12,6),sharey=True);labs=[]
  for ct in TYPES:
   for st in ['75-105d','330-400d']:
    row=d[d.celltype.eq(ct)&d.stratum.eq(st)].iloc[0];i=len(labs);labs.append(f'{ct} / {st} (n={int(row.n_pairs)})');col='#0072B2' if ct in ['B','MONO','NK'] else '#888888';axs[0].hlines(i,row.delta_CI_low,row.delta_CI_high,color=col);axs[0].plot(row.delta_ICC_e_minus_a,i,'o',color=col,ms=5);axs[1].hlines(i,row.ICC_repeat_p025,row.ICC_repeat_p975,color='#009E73');axs[1].plot(row.ICC_median,i,'o',color='#009E73');axs[1].plot(row.baseline_a_ICC,i,'d',color='#555555',ms=4,label='Baseline a' if i==0 else None)
  axs[0].set_yticks(np.arange(len(labs)),labs);axs[0].invert_yaxis();axs[0].axvline(0,color='black',ls='--',lw=.8);axs[0].set_xlabel('Paired donor-bootstrap 95% CI for ICC(e) - ICC(a)');axs[0].set_title('Primary e: full 68-donor internal extension\nBlue: prespecified F2 target types');axs[1].set_xlabel('ICC: median and 2.5–97.5 percentiles (50 repeats)');axs[1].set_title('Secondary c: repeated two-fold OOF\nSplit percentiles; not confidence intervals');axs[1].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.12))
  axs[1].set_xlim(0,1.02)
  for ax in axs:clean(ax);ax.grid(axis='x',alpha=.13)
  fig.tight_layout();save(fig,'Figure4_repair_confirmation')
 if (O/'Figure5_accuracy_bias_decomposition.csv').exists():
  a=pd.read_csv(O/'Figure5_accuracy_bias_decomposition.csv');fig=plt.figure(figsize=(15,10));grid=fig.add_gridspec(3,3,height_ratios=[1,1,1.1]);cols={'a':'#555555','c':'#009E73','e':'#D55E00'}
  for j,ds in enumerate(['OneK1K_test20','Allen_K0_full','AIDA_external']):
   for row,key,label in [(0,'MAE','MAE (years)'),(1,'Pearson_r','Pearson r')]:
    ax=fig.add_subplot(grid[row,j]);q=a[a.dataset.eq(ds)]
    for v,off in [('a',-.24),('c',0),('e',.24)]:vals=q[q.variant.eq(v)].set_index('celltype').reindex(TYPES)[key];ax.bar(np.arange(5)+off,vals,width=.22,color=cols[v],label=v)
    ax.set_xticks(np.arange(5),TYPES);ax.set_ylabel(label);clean(ax);ax.grid(axis='y',alpha=.12)
    if row==0:ax.set_title({'OneK1K_test20':'OneK1K holdout (20%)','Allen_K0_full':'Allen K0 (68 donors)','AIDA_external':'AIDA'}[ds])
    else:ax.set_ylim(-.2,1)
    if row==0 and j==0:legend_handles,legend_labels=ax.get_legend_handles_labels()
  ax=fig.add_subplot(grid[2,:]);q=a[a.dataset.eq('AIDA_external')].set_index(['celltype','variant']).reindex(pd.MultiIndex.from_product([TYPES,['a','c','e']]));fraction=q.bias_squared_fraction_MSE.to_numpy();ax.bar(np.arange(15),fraction,color='#CC79A7',label='Fixed offset squared / MSE');ax.bar(np.arange(15),1-fraction,bottom=fraction,color='#BBBBBB',label='Centered error variance / MSE');ax.set_ylim(0,1);ax.axhline(.5,color='black',ls=':',lw=.8);ax.set_xticks(np.arange(15),[ct+'/'+v for ct in TYPES for v in ['a','c','e']]);ax.set_ylabel('Fraction of AIDA mean squared error');ax.legend(frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,-.15));ax.set_title('AIDA error decomposition: diagnostic only; no target-outcome recalibration applied');clean(ax);fig.legend(legend_handles,legend_labels,frameon=False,title='Variant',ncol=3,loc='upper center',bbox_to_anchor=(.5,.995));fig.tight_layout(h_pad=2,rect=(0,0,1,.94));save(fig,'Figure5_accuracy_and_bias')
 print('ROUND4_PLOTS_COMPLETE',flush=True)
if __name__=='__main__':main()
