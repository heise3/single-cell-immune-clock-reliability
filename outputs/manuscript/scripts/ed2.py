import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
from matplotlib.lines import Line2D
rng=np.random.default_rng(20261025)
Dd=pd.read_csv(f'{DER}/Fig4a_within_person_change.csv')
SM=pd.read_csv(f'{DER}/Fig4b_pace_and_apparent_rejuvenation.csv')
fig=plt.figure(figsize=(183*MM,78*MM))
gs=fig.add_gridspec(1,2,width_ratios=[1.6,1],wspace=0.42)
fig.subplots_adjust(left=0.09,right=0.97,bottom=0.17,top=0.8)
ax=fig.add_subplot(gs[0])
for i,ct in enumerate(CT):
    for j,(coh,mk) in enumerate([('Allen','o'),('Terekhova','s')]):
        g=Dd[(Dd.cohort==coh)&(Dd.celltype==ct)]; xc=i+(j-0.5)*0.38; d=g.delta_pred.values
        jit=rng.uniform(-0.13,0.13,len(d)); neg=d<0
        ax.scatter(xc+jit[~neg],d[~neg],s=3.2,marker=mk,color=CTCOL[ct],alpha=0.35,lw=0)
        ax.scatter(xc+jit[neg],d[neg],s=3.2,marker=mk,color=CTCOL[ct],alpha=0.9,lw=0)
        ax.plot([xc-0.15,xc+0.15],[np.median(d)]*2,color=DARK,lw=0.9)
        ax.text(xc,-27.5,f'{int((d<=-2.75).sum())}',ha='center',fontsize=5.2,color=DARK)
ax.axhline(0,color=GREY,lw=0.5); ax.axhline(-2.75,color=RED,lw=0.5,ls=(0,(3,2)))
ax.text(1.01,-2.75,'−2.75 y',transform=ax.get_yaxis_transform(),fontsize=5.2,color=RED,va='center')
ax.text(-0.62,-27.5,'≤ −2.75 y:',fontsize=5.2,ha='right')
ax.set_xticks(range(5)); ax.set_xticklabels([CTLAB[c] for c in CT]); ax.set_xlim(-0.6,4.6); ax.set_ylim(-30,25)
ax.set_ylabel('Within-person change in\npredicted age (years)')
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=DARK,markersize=2.6,label='Allen, 330–400 d (n = 62)'),
                   Line2D([],[],marker='s',ls='',color=DARK,markersize=2.6,label='Terekhova, nominal 1–2 y (n = 101)')],
          loc='lower left',bbox_to_anchor=(0,1.02),ncol=2,fontsize=5.4,handletextpad=0.1)
axb=fig.add_subplot(gs[1])
for i,ct in enumerate(CT):
    for j,(coh,mk,fc) in enumerate([('Allen','o','white'),('Terekhova','s',None)]):
        r=SM[(SM.cohort==coh)&(SM.celltype==ct)].iloc[0]; y=i+(j-0.5)*0.3
        axb.plot([r.pace_CI_low,r.pace_CI_high],[y,y],color=CTCOL[ct],lw=0.8)
        axb.scatter(r.pace_years_per_year,y,s=12,marker=mk,color=fc or CTCOL[ct],edgecolor=CTCOL[ct],lw=0.8,zorder=3)
axb.axvline(1,color=DARK,lw=0.5,ls=(0,(3,2))); axb.axvline(0,color=GREY,lw=0.4)
axb.text(1.06,-0.85,'chronological-\ntime reference',fontsize=5,va='bottom')
axb.set_yticks(range(5)); axb.set_yticklabels([CTLAB[c] for c in CT]); axb.set_ylim(4.6,-1.0)
for t,c in zip(axb.get_yticklabels(),CT): t.set_color(CTCOL[c])
axb.set_xlim(-3.5,3); axb.set_xlabel('Mean change in predicted age per\nelapsed year (Allen) or nominal year (Terekhova)')
for lt,(axx,dx,dy) in zip('ab',[(ax,-0.13,0.09),(axb,-0.3,0.09)]):
    bb=axx.get_position(); fig.text(bb.x0+dx*bb.width,bb.y1+dy,lt,fontsize=8,fontweight='bold',va='bottom')
save(fig,'ED_Fig2_annual_followup')
