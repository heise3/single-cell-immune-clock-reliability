#!/usr/bin/env python3
"""Static manuscript plots from saved result CSVs; no inferential recalculation."""
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[2];O=P/'outputs/round5';CTS=['B','CD4T','CD8T','MONO','NK'];MODELS=['published','a','e'];COL={'published':'#738596','a':'#CF8743','e':'#277B78'};LABEL={'published':'Published','a':'Replica a','e':'Pseudobulk e'};S={'75-105d':'90 d','330-400d':'1 y'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
def save(fig,name):
 for ext in ['png','svg','pdf']:fig.savefig(O/(name+'.'+ext),dpi=220,bbox_inches='tight')
 plt.close(fig)
def main():
 z=pd.read_csv(O/'P1_three_model_ICC_long.csv');z=z[z.celltype.isin(['B','MONO','NK'])];z.to_csv(O/'FigureP1_three_models_plot_data.csv',index=False);delta=pd.read_csv(O/'P1_e_vs_published_paired_delta.csv');holm=pd.read_csv(O/'P1_F2_Holm6.csv');d=delta.merge(holm,on=['celltype','stratum']);d.to_csv(O/'FigureP1_differences_and_Holm_plot_data.csv',index=False);keys=[(ct,st) for ct in ['B','MONO','NK'] for st in S];fig,ax=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.2,1]});ys=np.arange(6)
 for i,(ct,st) in enumerate(keys):
  for j,m in enumerate(MODELS):
   r=z[z.celltype.eq(ct)&z.stratum.eq(st)&z.model.eq(m)].iloc[0];y=i+(j-1)*.19;ax[0].errorbar(r.ICC2_1,y,xerr=[[r.ICC2_1-r.ICC_CI_low],[r.ICC_CI_high-r.ICC2_1]],fmt='o',ms=4,color=COL[m],capsize=2,label=LABEL[m] if i==0 else None)
  r=d[d.celltype.eq(ct)&d.stratum.eq(st)].iloc[0];ax[1].errorbar(r.delta_e_minus_published,i,xerr=[[r.delta_e_minus_published-r.delta_CI_low],[r.delta_CI_high-r.delta_e_minus_published]],fmt='o',color=COL['e'],capsize=3);ax[1].text(.78,i,f'{r.Holm6_p:.3f} / {r.centered_Holm6_p_sensitivity:.3f}',va='center',fontsize=8)
 for a in ax:a.set_yticks(ys,[f'{ct} | {S[st]}' for ct,st in keys]);a.invert_yaxis();a.axvline(0,color='.6',lw=.7);a.grid(axis='x',alpha=.18)
 ax[0].set(xlabel='Age-acceleration ICC(2,1), pointwise 95% CI',xlim=(-.3,1.03),title='A  K0: published / replica a / e');ax[0].legend(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=3,fontsize=8,frameon=False);ax[1].set(xlabel='ICC(e) - ICC(published), paired 95% CI',xlim=(-.4,1.35),title='B  Published-model comparison');ax[1].text(.78,-.65,'Holm6 p vs a\ntail / centred sensitivity',fontsize=8,va='bottom');fig.tight_layout();save(fig,'FigureP1_model_comparison')
 r2=pd.read_csv(O/'P2_covariate_between_donor_R2.csv');age=pd.read_csv(O/'P2_chronological_age_signal.csv');ret=pd.read_csv(O/'P2_retained_reliability_gain.csv');fig,ax=plt.subplots(1,3,figsize=(13,4),gridspec_kw={'width_ratios':[1.1,.9,1.2]})
 for j,m in enumerate(MODELS):
  q=r2[r2.model.eq(m)].set_index('celltype').loc[CTS];ax[0].bar(np.arange(5)+(j-1)*.23,q.joint_R2,width=.22,color=COL[m],label=LABEL[m])
 ax[0].set(xticks=np.arange(5),xticklabels=CTS,ylabel='Between-donor R²: CMV + sex',ylim=(0,.46),title='A  Stable-covariate association');ax[0].legend(fontsize=7,frameon=False)
 mat=age.pivot(index='celltype',columns='model',values='Pearson_r').loc[CTS,MODELS];im=ax[1].imshow(mat,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto');ax[1].set(xticks=range(3),xticklabels=['Published','a','e'],yticks=range(5),yticklabels=CTS,title='B  Chronological-age Pearson r')
 for i in range(5):
  for j in range(3):ax[1].text(j,i,f'{mat.iloc[i,j]:.2f}',ha='center',va='center',color='white' if mat.iloc[i,j]>.65 else 'black')
 q=ret[ret.celltype.isin(['B','MONO','NK'])].set_index(['celltype','stratum']).loc[keys].reset_index()
 for i,r in q.iterrows():ax[2].plot([r.delta_e_minus_a_before,r.delta_e_minus_a_after],[i,i],color='.6');ax[2].scatter(r.delta_e_minus_a_before,i,color=COL['a'],s=25,label='Before adjustment' if i==0 else None);ax[2].scatter(r.delta_e_minus_a_after,i,color=COL['e'],s=25,label='CMV + sex adjusted' if i==0 else None)
 ax[2].set(yticks=range(6),yticklabels=[f'{c} | {S[s]}' for c,s in keys],xlabel='ICC(e) - ICC(a)',title='C  Gain after nested donor exclusion');ax[2].invert_yaxis();ax[2].axvline(0,lw=.7,color='.5');ax[2].legend(fontsize=7,frameon=False,loc='lower right');fig.tight_layout();save(fig,'FigureP2_covariate_validity');r2.to_csv(O/'FigureP2_R2_plot_data.csv',index=False);age.to_csv(O/'FigureP2_age_signal_plot_data.csv',index=False);q.to_csv(O/'FigureP2_retained_gain_plot_data.csv',index=False)
 vc=pd.read_csv(O/'P3_three_model_variance_components.csv');vc=vc[vc.metric.eq('age_adjusted')].copy();vc.to_csv(O/'FigureP3_variance_plot_data.csv',index=False);comps=['donor','visit_nested_plus_unmodelled_residual','batch','cell_sampling','pseudocell_algorithm'];colors=['#277B78','#ADB6C4','#B65D48','#D8AF61','#696099'];names=['Donor','Visit + unmodelled residual','Batch','Cell sampling','Pseudocell'];fig,ax=plt.subplots(1,2,figsize=(13,5),gridspec_kw={'width_ratios':[1.4,1]});positions=np.array([i*4+j for i in range(5) for j in range(3)]);q=vc.set_index(['celltype','model']).loc[[(c,m) for c in CTS for m in MODELS]];left=np.zeros(15)
 for c,col,n in zip(comps,colors,names):ax[0].barh(positions,q[c+'_fraction'],left=left,color=col,label=n);left+=q[c+'_fraction'].to_numpy()
 ax[0].set(yticks=positions,yticklabels=[f'{c} | {LABEL[m]}' for c in CTS for m in MODELS],xlabel='Fraction of total conditional variance',title='A  Age-acceleration variance fractions');ax[0].invert_yaxis();ax[0].legend(fontsize=7,loc='lower center',bbox_to_anchor=(.5,-.28),ncol=3,frameon=False)
 for j,m in enumerate(MODELS):
  t=vc[vc.model.eq(m)].set_index('celltype').loc[CTS];ax[1].bar(np.arange(5)+(j-1)*.23,t.batch_variance,width=.22,color=COL[m],label=LABEL[m])
 ax[1].set(xticks=range(5),xticklabels=CTS,ylabel='Conditional batch variance (years²)',title='B  Absolute batch variance');ax[1].legend(frameon=False);fig.tight_layout();save(fig,'FigureP3_repair_variance')
 f=O/'P4_external_technical_ICC.csv'
 if f.exists():
  z=pd.read_csv(f);fig,ax=plt.subplots(1,2,figsize=(11,4))
  for k,metric in enumerate(['prediction','age_acceleration']):
   for i,ct in enumerate(CTS):
    for j,m in enumerate(['published','e']):
     r=z[z.celltype.eq(ct)&z.clock.eq(m)&z.metric.eq(metric)].iloc[0];ax[k].errorbar(i+(j-.5)*.18,r.ICC2_1,yerr=[[r.ICC2_1-r.ICC_CI_low],[r.ICC_CI_high-r.ICC2_1]],fmt='o',color=COL[m],capsize=3,label=LABEL[m] if i==0 else None)
   ax[k].set(xticks=range(5),xticklabels=CTS,ylabel='ICC(2,1), donor bootstrap 95% CI',title=('Raw predicted age' if k==0 else 'LODO age acceleration'));ax[k].legend(frameon=False,fontsize=8);ax[k].axhline(0,color='.6',lw=.7)
  fig.suptitle('AIDA: 20 donors, same-draw technical GEX libraries (not longitudinal biology)',fontsize=10);fig.tight_layout();save(fig,'FigureP4_external_technical');z.to_csv(O/'FigureP4_external_technical_plot_data.csv',index=False)
if __name__=='__main__':main()
