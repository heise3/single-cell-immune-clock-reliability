#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[2];O=P/'outputs/round2_A1';r=pd.read_csv(O/'reliability_results.csv');s=r[r.clock.eq('sc-ImmuAging')];types=['B','CD4T','CD8T','MONO','NK'];strata=['7-9d','75-105d','330-400d'];colors=['#0072B2','#D55E00','#009E73'];fig,axs=plt.subplots(1,2,figsize=(10,5.6),sharey=True)
for ax,metric,title in zip(axs,['prediction','age_acceleration'],['Predicted age','LODO age acceleration']):
 for i,(stratum,color) in enumerate(zip(strata,colors)):
  d=s[(s.metric==metric)&(s.stratum==stratum)].set_index('celltype').reindex(types);y=np.arange(5)+(i-1)*.2;x=d.ICC2_1.to_numpy();lo=d.ICC_CI_low.to_numpy();hi=d.ICC_CI_high.to_numpy();ax.hlines(y,lo,hi,color=color,lw=1);ax.plot(x,y,'o',color=color,markersize=4,label=stratum)
 ax.axvline(.8,color='black',ls='--',lw=1);ax.axvline(0,color='#AAAAAA',lw=.5);ax.set_title(title);ax.set_xlabel('ICC(2,1), donor-bootstrap 95% CI');ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
axs[0].set_yticks(np.arange(5),types);axs[0].invert_yaxis();handles,labels=axs[1].get_legend_handles_labels();fig.legend(handles,labels,frameon=False,ncols=3,loc='lower center',bbox_to_anchor=(.5,0));fig.suptitle('Frozen K0 cohort: 68 donors / 320 steady-state visits',fontsize=12);fig.tight_layout(rect=(0,.065,1,1));fig.savefig(O/'scImmuAging_ICC_forest.svg');fig.savefig(O/'scImmuAging_ICC_forest.png',dpi=220);plt.close(fig)
v=pd.read_csv(O/'variance_decomposition.csv');v=v[(v.clock=='sc-ImmuAging')&(v.metric=='age_adjusted')].set_index('celltype').reindex(types);components=[('donor_fraction','Donor'),('batch_fraction','Batch'),('visit_nested_plus_unmodelled_residual_fraction','Visit + unmodelled residual'),('cell_sampling_fraction','Cell sampling'),('pseudocell_algorithm_fraction','Pseudocell randomness')];fig,ax=plt.subplots(figsize=(9,4));left=np.zeros(5)
for (key,label),color in zip(components,['#0072B2','#56B4E9','#999999','#E69F00','#CC79A7']):
 x=v[key].to_numpy();ax.barh(types,x,left=left,label=label,color=color);left+=x
assert np.allclose(left,1);ax.invert_yaxis();ax.set_xlabel('Fraction of conditional age-acceleration variance');ax.set_xlim(0,1);ax.legend(frameon=False,ncols=2,bbox_to_anchor=(0,-.2),loc='upper left');ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(O/'variance_decomposition_age_acceleration.svg');fig.savefig(O/'variance_decomposition_age_acceleration.png',dpi=220);plt.close(fig)
if (O/'RTM_null_results.csv').exists():
 rr=pd.read_csv(O/'RTM_null_results.csv');q=rr[rr.role=='PRIMARY'].iloc[0];a=np.load(O/('RTM_null_distributions_'+q.stratum+'.npz'));fig,axs=plt.subplots(1,2,figsize=(10,3.6))
 for ax,key,observed,label in [(axs[0],'negative_count',q.observed_n_negative,'Number of negative changes'),(axs[1],'negative_subgroup_median',q.observed_negative_median,'Median within negative-change subgroup (years)')]:
  z=a[key];z=z[np.isfinite(z)];bins=np.arange(z.min()-.5,z.max()+1.5,1) if key=='negative_count' else 35;ax.hist(z,bins=bins,color='#56B4E9',alpha=.8);ax.axvline(observed,color='#D55E00',lw=2,label='Observed');ax.set_xlabel(label);ax.set_ylabel('Null simulations');ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
 fig.suptitle('Prespecified null, CD8 age-acceleration SEM at approximately 90 days');fig.tight_layout();fig.savefig(O/'RTM_primary_null.svg');fig.savefig(O/'RTM_primary_null.png',dpi=220);plt.close(fig)
