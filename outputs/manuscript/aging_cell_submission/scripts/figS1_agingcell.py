# Existing Figure S1 source preserved; cosmetic legend placement only.
import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
from matplotlib.lines import Line2D
X=pd.read_csv(f'{DER}/ED1_prediction_vs_acceleration_ICC.csv')
mk={'7–9 d':'o','75–105 d':'D','330–400 d':'^','1–2 y (Terekhova)':'s'}
fig,ax=plt.subplots(figsize=(89*MM,116*MM))
fig.subplots_adjust(left=0.19,right=0.97,bottom=0.23,top=0.83)
ax.plot([-0.3,1],[-0.3,1],color=GREY,lw=0.5,ls=(0,(3,2)))
for _,r in X.iterrows():
    ax.scatter(r.ICC_prediction,r.ICC_age_acceleration,marker=mk[r.interval],s=16,color=CTCOL[r.celltype],edgecolor='white',lw=0.3,zorder=3)
ax.set_xlim(-0.3,1); ax.set_ylim(-0.3,1); ax.set_aspect('equal')
ax.set_xlabel('ICC of predicted age'); ax.set_ylabel('ICC of age acceleration')
colors=[Line2D([],[],marker='o',ls='',color=CTCOL[c],markersize=3.5,label=CTLAB[c]) for c in CT]
color_legend=ax.legend(handles=colors,loc='lower right',bbox_to_anchor=(1.0,1.035),
                       fontsize=5.4,ncol=2,handletextpad=0.2,columnspacing=0.8,borderaxespad=0)
ax.add_artist(color_legend)
shapes=[Line2D([],[],marker=m,ls='',color='#777777',markersize=3.5,label=k) for k,m in mk.items()]
interval_legend=ax.legend(handles=shapes,loc='upper center',bbox_to_anchor=(0.5,-0.20),
                         fontsize=5.4,ncol=2,handletextpad=0.2,columnspacing=0.8,borderaxespad=0)
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
assert not color_legend.get_window_extent(renderer).overlaps(ax.get_window_extent(renderer))
assert len(ax.collections)==len(X)==20
save(fig,'Figure_S1')
print('Figure S1: color legend upper right, outside axes; 20 points retained.')
