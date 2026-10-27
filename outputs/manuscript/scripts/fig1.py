import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
from matplotlib.patches import FancyBboxPatch, Circle
fig=plt.figure(figsize=(183*MM,116*MM))
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,183); ax.set_ylim(0,116); ax.axis('off')
def box(x,y,w,h,fc='white',ec=DARK,lw=0.5,r=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',fc=fc,ec=ec,lw=lw))
def txt(x,y,s,**k):
    k.setdefault('fontsize',6.2); k.setdefault('va','center'); k.setdefault('ha','left'); ax.text(x,y,s,**k)
ax.text(1,112,'a',fontsize=8,fontweight='bold',va='top')
txt(6,111.2,'Three layers of repeated measurement',fontsize=7,fontweight='bold',va='top')
cols=[(4,'Technical replicate','AIDA','20 donors · 36,805 cells · 21–66 y',
       ['Same blood draw, same PBMC pool','Two independent 10x libraries','Interval: 0 d']),
      (65,'Short-to-annual repeats','Allen Sound Life (GSE271896)','68 donors · 320 steady-state visits\n25–65 y',
       ['90 sequencing batches','Pairs: 7–9 d (33), 75–105 d (28),','330–400 d (62)']),
      (126,'External annual repeats','Terekhova blood atlas (HCA)','101 donors · 250 visits · 1.42M cells\n25–80 y',
       ['Earliest vs latest visit','Nominal 1 y (39) or 2 y (62)','97 of 101 pairs cross batch'])]
shades=['#F4F4F4','#ECECEC','#E3E3E3']
for (x,title,src,meta,lines),sh in zip(cols,shades):
    box(x,60,53,46,fc=sh,ec='none')
    txt(x+2.5,102,title,fontsize=6.8,fontweight='bold')
    txt(x+2.5,97.5,src,fontsize=6,color='#444444')
    txt(x+2.5,92,meta,fontsize=5.6,color='#444444',linespacing=1.5)
    # timeline glyph
    y0=79
    ax.plot([x+4,x+50],[y0,y0],color='#9A9A9A',lw=0.6)
    if x==4:
        for dx,c in [(12,DARK)]:
            ax.add_patch(Circle((x+dx,y0),1.4,fc=c,ec='none'))
        ax.plot([x+12,x+22],[y0,y0+3],color=DARK,lw=0.5); ax.plot([x+12,x+22],[y0,y0-3],color=DARK,lw=0.5)
        for yy in (y0+3,y0-3): box(x+22,yy-1.1,6,2.2,fc='white',ec=DARK,lw=0.5,r=0.4)
        txt(x+30,y0+3,'library 1',fontsize=5.4); txt(x+30,y0-3,'library 2',fontsize=5.4)
    elif x==65:
        pts=[(8,'#4D4D4D'),(10.5,'#4D4D4D'),(22,'#7A7A7A'),(42,'#B0B0B0')]
        for dx,c in pts: ax.add_patch(Circle((x+dx,y0),1.25,fc=c,ec='none'))
        for a,b,l,yo in [(8,10.5,'7–9 d',2.8),(10.5,22,'75–105 d',-3.6),(10.5,42,'330–400 d',3.6)]:
            ax.annotate('',xy=(x+b,y0+(1.6 if yo>0 else -1.6)),xytext=(x+a,y0+(1.6 if yo>0 else -1.6)),
                        arrowprops=dict(arrowstyle='-',color='#6E6E6E',lw=0.4,connectionstyle='arc3,rad=%s'%(-0.35 if yo>0 else 0.35)))
            txt(x+(a+b)/2,y0+yo+ (1.4 if yo>0 else -1.2),l,fontsize=5.2,ha='center',color='#4D4D4D')
    else:
        for dx in (10,26,42): ax.add_patch(Circle((x+dx,y0),1.25,fc='#6E6E6E',ec='none'))
        txt(x+18,y0+3.4,'year 1',fontsize=5.2,ha='center',color='#4D4D4D'); txt(x+34,y0+3.4,'year 2',fontsize=5.2,ha='center',color='#4D4D4D')
    for i,l in enumerate(lines): txt(x+2.5,70-i*3.6,l,fontsize=5.6)
# b pipeline
ax.text(1,54,'b',fontsize=8,fontweight='bold',va='top')
txt(6,53.2,'From cells to an individual-level reliability estimate',fontsize=7,fontweight='bold',va='top')
steps=[('PBMC\nscRNA-seq','one visit'),('Five lineages\nB · CD4$^+$ T · CD8$^+$ T\nmonocyte · NK','≥30 cells each'),
       ('sc-ImmuAging\nprediction','100 pseudocells\n× 15 cells'),('Age acceleration','prediction − f(age),\nleave-one-donor-out fit'),
       ('ICC(2,1) · SEM\nMDC$_{95}$','2,000 donor\nbootstraps')]
xs=[4,40,76,112,148]; w=27
for (s,sub),x in zip(steps,xs):
    box(x,23,w,22,fc='white',ec='#5A5A5A',lw=0.5)
    txt(x+w/2,37.5,s,fontsize=6,ha='center',fontweight='bold' if x in (112,148) else 'normal')
    txt(x+w/2,28,sub,fontsize=5.3,ha='center',color='#5A5A5A')
for a,b in zip(xs[:-1],xs[1:]):
    ax.annotate('',xy=(b-0.6,34),xytext=(a+w+0.6,34),arrowprops=dict(arrowstyle='-|>',color='#5A5A5A',lw=0.7,mutation_scale=8,shrinkA=0,shrinkB=0))
# c equation
txt(4,17,'Variance model',fontsize=6.2,fontweight='bold')
txt(4,11.5,r'$\mathrm{age\ acceleration}_{ijb}\;=\;\mu\;+\;u_i^{\,\mathrm{donor}}\;+\;v_b^{\,\mathrm{batch}}\;+\;e_{ij}^{\,\mathrm{visit}}\;+\;\varepsilon_{ij}^{\,\mathrm{cells}}\;+\;\varepsilon_{ij}^{\,\mathrm{pseudocell}}$',fontsize=6.6)
txt(4,5,r'$\mathrm{MDC}_{95}=1.96\,\sqrt{2}\,\mathrm{SEM}$: the smallest within-person change distinguishable from repeat-measurement variation with 95% confidence.',fontsize=5.8,color='#444444')
save(fig,'Fig1_design')
