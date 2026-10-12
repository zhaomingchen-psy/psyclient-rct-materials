import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
COL={'Lecture-only':'#7F7F7F','No-feedback':'#E69F00','Feedback':'#0072B2'}   # Okabe-Ito, grey control
ORDER=['Lecture-only','No-feedback','Feedback']
def setup():
    mpl.rcParams.update({
        'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans','DejaVu Sans'],
        'font.size':7,'axes.titlesize':8,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
        'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':2.5,'ytick.major.size':2.5,
        'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,
        'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':400,'figure.dpi':150})
def panel_label(ax,letter,x=-0.18,y=1.06):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom',ha='left')
