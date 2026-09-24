"""Minimal paired panels emphasizing out-of-range spouse correlations."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE.parents[1]/'figures'

def main():
    fits=json.loads((HERE/'theta1_unrestricted_vs_rm.json').read_text(encoding='utf-8'))['fits']
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(10,6.5),sharey=True,gridspec_kw={'wspace':.10})
    fig.subplots_adjust(left=.10,right=.95,top=.73,bottom=.18)
    accent='#aa4038'
    for ax,key,title,count in zip(axes,['unrestricted','constrained'],['Without constraint','With m = h²r'],[1,6]):
        values=sorted(f[key]['r'] for f in fits)
        ax.axhline(1,color=accent,lw=1.2,zorder=1)
        for i,r in enumerate(values):
            ax.scatter(i,r,s=65,color=accent if r>1 else '#a4abb0',edgecolors='white',linewidths=.7,zorder=3)
        ax.set_xlim(-.7,7.7);ax.set_ylim(0,4.5);ax.set_xticks([])
        ax.set_yticks([0,1,2,3,4]);ax.tick_params(axis='y',length=0,labelsize=11,pad=8)
        ax.spines[['top','right','bottom','left']].set_visible(False)
        ax.set_title(title,fontsize=14,pad=17)
        ax.text(.5,-.10,f'{count} of 8 above 1',transform=ax.transAxes,ha='center',fontsize=16,color=accent,fontweight='bold')
        ax.text(.5,-.17,'Each dot is one trait; sorted within each panel',transform=ax.transAxes,ha='center',fontsize=8,color='#747d84')
    axes[0].set_ylabel('Estimated spouse correlation r',fontsize=12,labelpad=13)
    axes[1].tick_params(labelleft=False)
    axes[1].text(7.6,1.10,'Maximum possible: 1',ha='right',va='bottom',fontsize=10,color=accent)
    axes[1].annotate('Company director: 4.09',(7,max(f['constrained']['r'] for f in fits)),
                     xytext=(-10,15),textcoords='offset points',ha='right',fontsize=10,color=accent)
    fig.text(.10,.94,'One impossible correlation becomes six.',fontsize=23,fontweight='bold',ha='left',color='#20252a')
    fig.text(.10,.865,'Imposing the phenotypic-assortment constraint barely moves heritability\nbut pushes most spouse correlations above their mathematical limit.',fontsize=12,ha='left',color='#545e66',linespacing=1.5)
    fig.text(.10,.035,'Table 4.5 · Both θ = 1 · Equal-weight log residuals · No upper bounds imposed',fontsize=9,color='#747d84')
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'book_theta1_spouse_r_minimal.{ext}',dpi=180,facecolor='white')
    plt.close(fig)

if __name__=='__main__':main()
