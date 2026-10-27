import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt, os, numpy as np, pandas as pd
from pathlib import Path
from matplotlib import font_manager
_font_root=Path(os.environ.get('IMMUNE_CLOCK_FONT_DIR', '/usr/share/fonts/truetype/liberation2'))
for _font in _font_root.glob('LiberationSans-*.ttf'):
    font_manager.fontManager.addfont(str(_font))
mpl.rcParams.update({
 'font.family':'Liberation Sans','font.size':6.5,'axes.titlesize':7,'axes.labelsize':6.5,
 'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,'legend.frameon':False,
 'axes.linewidth':0.5,'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2.5,'ytick.major.size':2.5,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none',
 'lines.linewidth':0.8,'patch.linewidth':0.5,'axes.titlepad':4,'figure.dpi':150,'savefig.dpi':450,
 'mathtext.default':'regular'})
_HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.abspath(os.path.join(_HERE,'..','..'))
OUT=os.path.abspath(os.path.join(_HERE,'..'))
FIG=OUT+'/figures'; DER=OUT+'/source_data'
CT=['B','CD4T','CD8T','MONO','NK']
CTLAB={'B':'B','CD4T':'CD4$^+$ T','CD8T':'CD8$^+$ T','MONO':'Monocyte','NK':'NK'}
CTCOL={'B':'#0072B2','CD4T':'#E69F00','CD8T':'#D55E00','MONO':'#CC79A7','NK':'#009E73'}
GREY='#8C8C8C'; DARK='#2B2B2B'; LIGHT='#EDEDED'; RED='#B2182B'
MM=1/25.4
def panel(ax,letter,x=-0.12,y=1.04):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=8,fontweight='bold',va='bottom',ha='left')
def save(fig,name):
    # Layout QA only. Cached scientific source CSVs are never written here.
    from pathlib import Path
    import json
    from matplotlib.text import Text
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    canvas=fig.bbox
    outside=[]
    for text in fig.findobj(match=Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        bb=text.get_window_extent(renderer)
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > canvas.x1+0.5 or bb.y1 > canvas.y1+0.5:
            outside.append(text.get_text())
    qa=Path(OUT)/'layout_QA'
    qa.mkdir(exist_ok=True)
    (qa/(name+'.json')).write_text(json.dumps({'figure':name,'width_mm':fig.get_figwidth()*25.4,'height_mm':fig.get_figheight()*25.4,'text_outside_fixed_canvas':outside,'layout_only':True},ensure_ascii=False,indent=2))
    for ext in ['pdf','svg','png']:
        fig.savefig(f'{FIG}/{name}.{ext}',bbox_inches=None)
    plt.close(fig)
