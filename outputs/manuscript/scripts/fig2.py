import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
L=pd.read_csv(f'{DER}/Fig2a_ICC_by_condition.csv')
M=pd.read_csv(f'{DER}/Fig2b_MDC95.csv')
V=pd.read_csv(f'{DER}/Fig2c_age_validity.csv')
conds=['Technical\nreplicate','7–9 d\nsame batch','7–9 d\ncross batch','75–105 d','330–400 d','1–2 y\n(external)']
ns=['n = 20','n = 18','n = 15','n = 28','n = 62','n = 101']

fig=plt.figure(figsize=(183*MM,148*MM))
gs=fig.add_gridspec(2,2,height_ratios=[1.05,1.0],width_ratios=[1.25,1],hspace=0.90,wspace=0.82)
fig.subplots_adjust(left=0.105,right=0.91,bottom=0.12,top=0.86)
# ---- a: dot matrix
ax=fig.add_subplot(gs[0,:])
cmap=LinearSegmentedColormap.from_list('icc',['#F6EFE6','#E3CBB0','#9DB7C4','#4A7C96','#163E57'])
norm=Normalize(0,1)
for _,r in L.iterrows():
    i=CT.index(r.celltype); k=r.condition_index; v=r.ICC
    size=55+215*max(v,0)
    ax.scatter(k,i,s=size,color=cmap(norm(max(v,0))),edgecolor='white',lw=0.6,zorder=3)
    ax.text(k,i,f'{v:.2f}',ha='center',va='center',fontsize=5.8,color='white' if v>0.55 else DARK,zorder=4)
ax.set_xlim(-0.6,5.6); ax.set_ylim(4.6,-0.9)
ax.set_yticks(range(5)); ax.set_yticklabels([CTLAB[c] for c in CT])
for t,c in zip(ax.get_yticklabels(),CT): t.set_color(CTCOL[c]); t.set_fontweight('bold')
ax.set_xticks(range(6)); ax.set_xticklabels(conds,fontsize=6)
ax.tick_params(axis='x',length=0,pad=3); ax.tick_params(axis='y',length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.axvline(1.5,color='#BDBDBD',lw=0.5,ls=(0,(2,2)),ymin=0.02,ymax=0.86)
for k,n in enumerate(ns): ax.text(k,-0.72,n,ha='center',va='center',fontsize=5.3,color='#6A6A6A')
ax.annotate('',xy=(1.42,-1.25),xytext=(-0.42,-1.25),xycoords='data',arrowprops=dict(arrowstyle='-',color=DARK,lw=0.5),annotation_clip=False)
ax.annotate('',xy=(5.42,-1.25),xytext=(1.58,-1.25),xycoords='data',arrowprops=dict(arrowstyle='-',color=DARK,lw=0.5),annotation_clip=False)
ax.text(0.5,-1.42,'Within a batch',ha='center',va='bottom',fontsize=6.3,fontweight='bold')
ax.text(3.5,-1.42,'Across batches and visits',ha='center',va='bottom',fontsize=6.3,fontweight='bold')
ax.set_ylim(4.6,-1.1)
cax=fig.add_axes([0.920,0.61,0.009,0.22]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label('ICC(2,1),\nage acceleration',fontsize=5.6); cb.ax.tick_params(labelsize=5.3,width=0.4,length=2); cb.outline.set_linewidth(0.4)
panel(ax,'a',x=-0.07,y=1.1)
# ---- b: MDC
axb=fig.add_subplot(gs[1,0])
cl=[('Technical replicate',None),('7–9 d','7-9d'),('75–105 d','75-105d'),('330–400 d','330-400d'),('1–2 y (external)',None)]
mk=['D','o','o','o','s']; edge=['#B5B5B5','#8A8A8A','#5C5C5C','#2E2E2E','#000000']
w=0.15
axb.axhspan(0,2.75,color='#F7E4E1',lw=0,zorder=0)
axb.axhline(2.75,color=RED,lw=0.5,ls=(0,(3,2)),zorder=1)
axb.text(-0.45,0.65,'Post-BCG CD8$^+$ responder median: 2.75 y',fontsize=5.1,color=RED,va='bottom')
for j,(lab,_) in enumerate(cl):
    d=M[M.condition==lab]
    for i,ct in enumerate(CT):
        r=d[d.celltype==ct].iloc[0]; x=i+(j-2)*w
        axb.plot([x,x],[r.CI_low,r.CI_high],color=CTCOL[ct],lw=0.75,zorder=2)
        axb.scatter(x,r.MDC95,marker=mk[j],s=11,color=CTCOL[ct],edgecolor=edge[j],lw=0.55,zorder=3)
axb.set_xticks(range(5)); axb.set_xticklabels([CTLAB[c] for c in CT])
axb.set_ylabel('MDC$_{95}$ of age acceleration (years)'); axb.set_ylim(0,22); axb.set_xlim(-0.55,4.55)
h=[Line2D([],[],marker=mk[j],ls='',color='#D9D9D9',markeredgecolor=edge[j],markersize=3.4,label=cl[j][0]) for j in range(5)]
axb.legend(handles=h,ncol=2,loc='lower left',bbox_to_anchor=(0,1.10),handletextpad=0.1,columnspacing=0.6,fontsize=5.5)
panel(axb,'b',x=-0.19,y=1.26)
# ---- c: validity
axc=fig.add_subplot(gs[1,1])
axc.axvspan(-0.3,0.5,color='#F4F4F4',lw=0,zorder=0)
for i,ct in enumerate(CT):
    a=V[(V.celltype==ct)&(V.cohort=='Allen')].iloc[0]; b=V[(V.celltype==ct)&(V.cohort=='Terekhova')].iloc[0]
    axc.plot([b.CI_low,b.CI_high],[i+0.13,i+0.13],color=CTCOL[ct],lw=0.8)
    axc.scatter(b.pearson_r,i+0.13,s=13,marker='s',color=CTCOL[ct],zorder=3)
    axc.scatter(a.pearson_r,i-0.13,s=13,marker='o',color='white',edgecolor=CTCOL[ct],lw=0.8,zorder=3)
axc.axvline(0,color=GREY,lw=0.4)
axc.set_yticks(range(5)); axc.set_yticklabels([CTLAB[c] for c in CT]); axc.set_ylim(4.6,-0.6)
for t,c in zip(axc.get_yticklabels(),CT): t.set_color(CTCOL[c])
axc.set_xlim(-0.3,1); axc.set_xlabel('Donor-level r, predicted vs chronological age')
h=[Line2D([],[],marker='o',ls='',color='white',markeredgecolor=DARK,markersize=3.4,label='Allen (n = 68)'),Line2D([],[],marker='s',ls='',color=DARK,markersize=3.4,label='Terekhova (n = 101)')]
axc.legend(handles=h,ncol=1,loc='lower left',bbox_to_anchor=(0,1.10),handletextpad=0.1,columnspacing=0.8,fontsize=5.5)
panel(axc,'c',x=-0.34,y=1.26)
save(fig,'Fig2_reliability')
print(V.round(3).to_string())
