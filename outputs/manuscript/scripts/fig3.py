import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
from matplotlib.lines import Line2D
ram_source=pd.read_csv(f'{DER}/Fig3a_RAM_same_cross.csv')
pub_source=pd.read_csv(f'{DER}/Fig2a_ICC_by_condition.csv')
pub_source=pub_source[pub_source.condition_index.isin([1,2])].copy()
pub_source['clock']='sc-ImmuAging'
pub_source['batch_relation']=pub_source.condition_index.map({1:'same_batch',2:'cross_batch'})
pub_source['ICC2_1']=pub_source.ICC
s7=pd.concat([ram_source,pub_source],ignore_index=True)
vc=pd.read_csv(f'{DER}/Fig3b_variance_components.csv').set_index('celltype')
ba=pd.read_csv(f'{DER}/Fig3c_batch_adjustment.csv')
fig=plt.figure(figsize=(183*MM,102*MM))
gs=fig.add_gridspec(1,3,width_ratios=[0.9,1.15,1.2],wspace=0.95)
fig.subplots_adjust(left=0.095,right=0.88,bottom=0.22,top=0.74)
# a
ax=fig.add_subplot(gs[0])
ram=s7[s7.clock.str.startswith('RAM')]
for (clk,ctt),g in ram.groupby(['clock','celltype']):
    g=g.set_index('batch_relation')
    ax.plot([0,1],[g.loc['same_batch','ICC2_1'],g.loc['cross_batch','ICC2_1']],color='#C8C8C8',lw=0.45,zorder=1)
for ct in CT:
    g=s7[(s7.clock=='sc-ImmuAging')&(s7.celltype==ct)].set_index('batch_relation')
    y=[g.loc['same_batch','ICC2_1'],g.loc['cross_batch','ICC2_1']]
    ax.plot([0,1],y,color=CTCOL[ct],lw=0.9,zorder=3); ax.scatter([0,1],y,s=12,color=CTCOL[ct],zorder=4)
    lab_y={'B':0.47,'NK':0.415,'MONO':0.36}.get(ct,y[1])
    ax.text(1.08,lab_y,CTLAB[ct],color=CTCOL[ct],fontsize=5.6,va='center')
ax.set_xticks([0,1]); ax.set_xticklabels(['Same\nbatch','Different\nbatch']); ax.set_xlim(-0.25,1.55)
ax.set_ylim(0,1.0); ax.set_ylabel('ICC, 7–9-day repeat')
ax.text(0.5,-0.24,'n = 18 vs 15 donor pairs',transform=ax.transAxes,ha='center',fontsize=5.3,color='#6A6A6A')
ax.legend(handles=[Line2D([],[],color='#C8C8C8',lw=0.6,label='RAM signatures (16)')],bbox_to_anchor=(-0.05,1.04),loc='lower left',fontsize=5.3,handlelength=1.2)

# Cached figure source data remain unchanged.
# b
axb=fig.add_subplot(gs[1])
comps=[('donor_fraction','Stable between-donor','#1F4E67'),('batch_fraction','Sequencing batch','#D98C3F'),
       ('visit_nested_plus_unmodelled_residual_fraction','Visit and residual','#B9B9B9'),
       ('cell_sampling_fraction','Cell sampling','#7FA9BE'),('pseudocell_algorithm_fraction','Pseudocell draw','#E6D3A3')]
for i,ct in enumerate(CT):
    left=0
    for col,lab,c in comps:
        v=vc.loc[ct,col]; axb.barh(i,v,left=left,color=c,height=0.62,edgecolor='white',lw=0.4,label=lab if i==0 else None)
        if v>0.08: axb.text(left+v/2,i,f'{v*100:.0f}',ha='center',va='center',fontsize=5.3,color='white' if col in ('donor_fraction',) else DARK)
        left+=v
axb.set_yticks(range(5)); axb.set_yticklabels([CTLAB[c] for c in CT]); axb.set_ylim(4.5,-0.5)
for t,c in zip(axb.get_yticklabels(),CT): t.set_color(CTCOL[c])
axb.set_xlim(0,1); axb.set_xticks([0,0.25,0.5,0.75,1]); axb.set_xticklabels(['0','25','50','75','100'])
axb.set_xlabel('Share of variance in\nage acceleration (%)')
axb.legend(ncol=1,loc='lower left',bbox_to_anchor=(-0.02,1.04),fontsize=5.3,handlelength=0.9,handletextpad=0.3,columnspacing=0.6)
axb.tick_params(axis='y',length=0)

# Cached figure source data remain unchanged.
# c
axc=fig.add_subplot(gs[2])
strata=[('7-9d','7–9 d'),('75-105d','75–105 d'),('330-400d','330–400 d')]
ypos=[];labs=[]
y=0
for st,sl in strata:
    for ct in CT:
        r=ba[(ba.celltype==ct)&(ba.stratum==st)].iloc[0]
        axc.annotate('',xy=(r.ICC_batch_adjusted,y),xytext=(r.ICC_before,y),arrowprops=dict(arrowstyle='-|>',color=CTCOL[ct],lw=0.7,mutation_scale=5,shrinkA=0,shrinkB=1.5))
        axc.scatter(r.ICC_before,y,s=10,color='white',edgecolor=CTCOL[ct],lw=0.7,zorder=3)
        y+=1
    axc.text(1.02,y-3,sl,fontsize=5.8,va='center',ha='left',transform=axc.get_yaxis_transform())
    y+=0.8
yt=[];yl=[];y=0
for st,sl in strata:
    for ct in CT: yt.append(y); yl.append(CTLAB[ct]); y+=1
    y+=0.8
axc.set_yticks(yt); axc.set_yticklabels(yl,fontsize=5.2); axc.set_ylim(y-0.6,-0.8); axc.tick_params(axis='y',length=0)
for t,lab in zip(axc.get_yticklabels(),yl):
    ct=[c for c in CT if CTLAB[c]==lab][0]; t.set_color(CTCOL[ct])
axc.set_xlim(0,1); axc.set_xlabel('ICC of age acceleration')
for yy in [4.9,10.7]: axc.axhline(yy,color='#E0E0E0',lw=0.5)
axc.legend(handles=[Line2D([],[],marker='o',ls='',color='white',markeredgecolor=DARK,markersize=3,label='observed'),
                    Line2D([],[],marker=r'$\rightarrow$',ls='',color=DARK,markersize=5,label='after removing estimated batch effect')],
           loc='lower left',bbox_to_anchor=(-0.05,1.0),fontsize=5.3,ncol=1,handletextpad=0.3)

# Cached figure source data remain unchanged.
for lt,(axx,dx) in zip('abc',[(ax,-0.42),(axb,-0.25),(axc,-0.3)]):
    bb=axx.get_position(); fig.text(bb.x0+dx*bb.width,0.99,lt,fontsize=8,fontweight='bold',va='top')
save(fig,'Fig3_batch')
