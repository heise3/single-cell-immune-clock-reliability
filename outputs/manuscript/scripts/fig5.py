import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
from matplotlib.lines import Line2D
p1=pd.read_csv(f'{DER}/Fig5ab_three_model_ICC.csv')
p2=pd.read_csv(f'{DER}/Fig5d_age_validity.csv')
p4=pd.read_csv(f'{DER}/Fig5c_technical_ICC.csv'); p4=p4[p4.metric=='age_acceleration']
M={'published':('Published sc-ImmuAging','o','white',DARK),'a':('Matched retrained (a)','^','#BDBDBD','#7A7A7A'),'e':('Donor pseudobulk (e)','s',None,None)}
fig=plt.figure(figsize=(183*MM,147*MM))
gs=fig.add_gridspec(2,2,wspace=0.67,hspace=0.83)
fig.subplots_adjust(left=0.14,right=0.95,bottom=0.12,top=0.85)
def dots(ax,df,col,lo,hi,models,off=0.22):
    for i,ct in enumerate(CT):
        for j,m in enumerate(models):
            r=df[(df.celltype==ct)&(df.model==m)]
            if r.empty: continue
            r=r.iloc[0]; y=i+(j-(len(models)-1)/2)*off
            lab,mk,fc,ec=M[m]; fc=fc or CTCOL[ct]; ec=ec or CTCOL[ct]
            if lo and np.isfinite(r[lo]): ax.plot([r[lo],r[hi]],[y,y],color=CTCOL[ct] if m=='e' else '#9A9A9A',lw=0.7)
            ax.scatter(r[col],y,marker=mk,s=12,facecolor=fc,edgecolor=ec,lw=0.7,zorder=3)
    ax.set_yticks(range(5)); ax.set_ylim(4.6,-0.6); ax.tick_params(axis='y',length=0)
for k,(st,title) in enumerate([('75-105d','Biological repeat\n75–105 d'),('330-400d','Biological repeat\n330–400 d')]):
    ax=fig.add_subplot(gs[k]); dots(ax,p1[p1.stratum==st],'ICC2_1','ICC_CI_low','ICC_CI_high',['published','a','e'])
    ax.set_xlim(-0.3,1); ax.set_title(title,fontsize=6.2); ax.set_xlabel('ICC of age acceleration')
    ax.set_yticklabels([CTLAB[c] for c in CT] if k==0 else [])
    if k==0:
        for t,c in zip(ax.get_yticklabels(),CT): t.set_color(CTCOL[c])
    ax.axvline(0,color=GREY,lw=0.4)
ax3=fig.add_subplot(gs[2]); dots(ax3,p4.rename(columns={'clock':'model'}),'ICC2_1','ICC_CI_low','ICC_CI_high',['published','e'],off=0.25)
ax3.set_xlim(0.3,1); ax3.set_title('Technical replicate (AIDA)',fontsize=6.2); ax3.set_xlabel('ICC of age acceleration'); ax3.set_yticklabels([CTLAB[c] for c in CT])
for t,c in zip(ax3.get_yticklabels(),CT): t.set_color(CTCOL[c])
ax4=fig.add_subplot(gs[3]); dots(ax4,p2,'Pearson_r',None,None,['published','a','e'])
ax4.set_xlim(-0.1,1); ax4.set_title('Age validity (Allen, n = 68)',fontsize=6.2); ax4.set_xlabel('Donor-level r with\nchronological age'); ax4.set_yticklabels([])
ax4.axvline(0,color=GREY,lw=0.4)
h=[Line2D([],[],marker=M[m][1],ls='',markerfacecolor=M[m][2] or '#555555',markeredgecolor=M[m][3] or '#555555',markersize=3.6,label=M[m][0]) for m in ['published','a','e']]
fig.legend(handles=h,ncol=3,loc='upper center',bbox_to_anchor=(0.54,0.98),fontsize=5.8,handletextpad=0.2,columnspacing=1.4)
axs=[fig.axes[0],fig.axes[1],ax3,ax4]
for lt,axx in zip('abcd',axs):
    bb=axx.get_position(); fig.text(bb.x0-0.08,bb.y1+0.055,lt,fontsize=8,fontweight='bold',va='bottom')
save(fig,'Fig5_aggregation_tradeoff')
