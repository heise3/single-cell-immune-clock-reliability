import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
from matplotlib.lines import Line2D
rng=np.random.default_rng(20261025)
M=pd.read_csv(f'{DER}/Fig4a_matched_interval_change.csv')
S=pd.read_csv(f'{DER}/Fig4a_matched_interval_summary.csv')
bcg=pd.read_csv(f'{DER}/Fig4c_BCG_38_paired_CD8.csv')
curve=pd.read_csv(f'{DER}/Fig4d_SEM_sensitivity_curve.csv')
crit=pd.read_csv(f'{DER}/Fig4d_critical_SEM_vs_CD8_estimates.csv'); crit=crit[crit.metric=='age_acceleration']
fig=plt.figure(figsize=(183*MM,140*MM))
gs=fig.add_gridspec(2,3,height_ratios=[1,1],hspace=0.78,wspace=0.6)
fig.subplots_adjust(left=0.09,right=0.95,bottom=0.2,top=0.9)
# a ---------------------------------------------------------------
ax=fig.add_subplot(gs[0,:])
A=M[(M.source=='Allen')&(M.stratum=='75-105d')]
xs={ct:i for i,ct in enumerate(CT)}; xb=6.0; lab_n={}
for ct in CT:
    g=A[A.celltype==ct]; d=g.change.values; j=rng.uniform(-0.17,0.17,len(d)); neg=d<0
    ax.scatter(xs[ct]+j[~neg],d[~neg],s=6,color=CTCOL[ct],alpha=0.35,lw=0)
    ax.scatter(xs[ct]+j[neg],d[neg],s=6,color=CTCOL[ct],alpha=0.95,lw=0)
    ax.plot([xs[ct]-0.22,xs[ct]+0.22],[np.median(d)]*2,color=DARK,lw=0.9)
    s=S[(S.source=='Allen')&(S.stratum=='75-105d')&(S.celltype==ct)].iloc[0]
    lab_n[ct]=f'{int(s.n_decline)}/{int(s.n)}'
d=bcg.change.values; j=rng.uniform(-0.17,0.17,len(d)); neg=d<0
ax.scatter(xb+j[~neg],d[~neg],s=7,marker='s',color='#BDBDBD',lw=0)
ax.scatter(xb+j[neg],d[neg],s=7,marker='s',color=CTCOL['CD8T'],lw=0)
ax.plot([xb-0.22,xb+0.22],[np.median(d)]*2,color=DARK,lw=0.9)
s=S[S.source=='BCG'].iloc[0]; lab_n['BCG']=f'{int(s.n_decline)}/{int(s.n)}'
ax.axvline(5.0,color='#CFCFCF',lw=0.5)
ax.axhline(0,color=GREY,lw=0.5); ax.axhline(-2.75,color=RED,lw=0.5,ls=(0,(3,2)))
ax.text(5.08,-3.2,'−2.75 y',fontsize=5.2,color=RED,va='top',ha='left')
ax.set_xticks(list(range(5))+[xb]); ax.set_xticklabels([f'{CTLAB[c]}\n{lab_n[c]} declined' for c in CT]+[f'CD8$^+$ T\n{lab_n["BCG"]} declined'])
ax.set_xlim(-0.6,6.5); ax.set_ylim(-20,24.5); ax.set_ylabel('Within-person change in\npredicted age (years)')
ax.text(2.0,26.5,'Healthy adults, 78–105 days apart (Allen, n = 28)\nNo recorded vaccination between draws',ha='center',fontsize=6,fontweight='bold')
ax.text(xb,26.5,'3 months after BCG\n(Li et al., 2025; n = 38)',ha='center',fontsize=6,fontweight='bold')
# b ---------------------------------------------------------------
axb=fig.add_subplot(gs[1,0])
bins=np.arange(-12,26,2)
axb.hist(bcg.change,bins=bins,color='#D9D9D9',edgecolor='white',lw=0.4)
axb.hist(bcg.change[bcg.change<0],bins=bins,color=CTCOL['CD8T'],edgecolor='white',lw=0.4)
xx=np.linspace(-12,24,300); mu=bcg.change.mean()
for sem,ls,lab in [(3.377,'-','SEM 3.38 y (75–105 d)'),(2.081,(0,(2,1.5)),'SEM 2.08 y (7–9 d)')]:
    sd=np.sqrt(2)*sem; pdf=np.exp(-0.5*((xx-mu)/sd)**2)/(sd*np.sqrt(2*np.pi))
    axb.plot(xx,pdf*len(bcg)*2,color=DARK,lw=0.7,ls=ls,label=lab)
axb.axvline(0,color=GREY,lw=0.4)
axb.set_xlabel('Change in CD8$^+$ T predicted age\nafter BCG (years)'); axb.set_ylabel('Donors')
axb.legend(loc='upper left',bbox_to_anchor=(0,-0.36),fontsize=5.1,handlelength=1.6,title='Measurement-error model',title_fontsize=5.1)
# c ---------------------------------------------------------------
axc=fig.add_subplot(gs[1,1])
axc.plot(curve.SEM,curve.negative_count_Holm_p,color='#5A5A5A',lw=0.8,label='number declining')
axc.plot(curve.SEM,curve.negative_median_Holm_p,color='#5A5A5A',lw=0.8,ls=(0,(2,1.5)),label='median decline')
axc.axhline(0.05,color=RED,lw=0.5); axc.text(5.9,0.07,'P = 0.05',fontsize=5.1,color=RED,ha='right',va='bottom')
axc.axvline(1.8,color=RED,lw=0.5,ls=(0,(3,2))); axc.set_title('Critical SEM: 1.8 y',fontsize=5.8,color=RED,pad=6)
for _,r in crit.iterrows():
    lab={'7-9d':'7–9 d','75-105d':'75–105 d','330-400d':'330–400 d'}[r.stratum]
    yv={'7-9d':0.86,'75-105d':0.72,'330-400d':0.58}[r.stratum]
    axc.plot([r.SEM_CI_low,r.SEM_CI_high],[yv,yv],color=CTCOL['CD8T'],lw=0.9); axc.scatter(r.SEM,yv,s=10,color=CTCOL['CD8T'],zorder=3)
    axc.text(1.03,yv,lab,transform=axc.get_yaxis_transform(),fontsize=5.1,va='center')
axc.set_xlim(0.5,6); axc.set_ylim(0,1.02); axc.set_xlabel('Assumed repeat-measurement\nSEM (years)'); axc.set_ylabel('Holm-adjusted P')
axc.legend(loc='upper left',bbox_to_anchor=(0,-0.36),fontsize=5.1,handlelength=1.6)
# d ---------------------------------------------------------------
axd=fig.add_subplot(gs[1,2])
axd.scatter(bcg.baseline,bcg.change,s=7,color=np.where(bcg.change<0,CTCOL['CD8T'],'#BDBDBD'),lw=0)
b,a=np.polyfit(bcg.baseline,bcg.change,1); xx=np.linspace(bcg.baseline.min(),bcg.baseline.max(),10); axd.plot(xx,a+b*xx,color=DARK,lw=0.7)
r=np.corrcoef(bcg.baseline,bcg.change)[0,1]
axd.set_title(f'r = {r:.2f}'.replace('-','−'),fontsize=5.8,pad=6)
axd.axhline(0,color=GREY,lw=0.4)
axd.set_xlabel('Baseline CD8$^+$ T\npredicted age (years)'); axd.set_ylabel('Change after BCG (years)')
for lt,(axx,dx,dy) in zip('abcd',[(ax,-0.075,0.06),(axb,-0.33,0.05),(axc,-0.33,0.05),(axd,-0.33,0.05)]):
    bb=axx.get_position(); fig.text(bb.x0+dx*bb.width,bb.y1+dy,lt,fontsize=8,fontweight='bold',va='bottom')
save(fig,'Figure_4')
