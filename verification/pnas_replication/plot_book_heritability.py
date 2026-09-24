"""Compare exp(a) with PH constraint-based h2 for book Table 4.5 OLS.

Maintained reporting hypothesis, not a reconstruction of Clark's implementation.
Run book_table45_ols.py first. No estimates clipped or confidence intervals claimed.
"""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE.parents[1]/'figures'

def main():
    report=json.loads((HERE/'book_table45_ols_results.json').read_text(encoding='utf-8'))
    fits=report['fits'];OUT.mkdir(exist_ok=True)
    labels=['House value','IMD','Company director','Occupation 1780','Occupation 1860',
            'Education 1780','Education 1860','Literacy']
    # offsets in points are explicit so closely spaced labels do not obscure points.
    offsets=[[(5,29),(-15,-20),(12,-18),(-12,-20),(30,0),(8,-20),(10,-25),(8,10)],
             [(-12,-24),(-10,-22),(10,9),(-12,-20),(10,4),(10,-18),(-10,-23),(8,9)]]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(13,6.8),sharex=True,sharey=True)
    fig.subplots_adjust(left=.075,right=.985,top=.79,bottom=.25,wspace=.12)
    titles=['A. Treating exp(a) as heritability','B. Recovering heritability using m = h²r']
    for panel,(ax,key,title) in enumerate(zip(axes,['theta_h2','h2'],titles)):
        ax.axvspan(1,2.3,color='#f9eded',zorder=0)
        ax.axvline(1,color='#b65252',lw=1.1,ls='--',zorder=1)
        ax.set_xlim(0,2.3);ax.set_ylim(.68,.94)
        ax.set_xticks(np.arange(0,2.26,.25));ax.set_yticks([.70,.75,.80,.85,.90])
        ax.grid(axis='y',color='#e0e5ea',lw=.7)
        ax.set_title(title,loc='left',fontsize=12,pad=15,fontweight='bold')
        for i,(f,label) in enumerate(zip(fits,labels)):
            invalid=bool(f['bound_violations'])
            ax.scatter(f[key],f['b'],s=65,c='#b7463e' if invalid else '#286c91',
                       marker='X' if invalid else 'o',edgecolors='white',linewidths=.7,zorder=4)
            dx,dy=offsets[panel][i]
            ax.annotate(label,(f[key],f['b']),xytext=(dx,dy),textcoords='offset points',
                        ha='right' if dx<0 else 'left',va='center',fontsize=9,
                        arrowprops=dict(arrowstyle='-',color='#a5adb5',lw=.65),zorder=5)
        ax.set_xlabel('exp(a) = θ × h²' if panel==0 else 'h² = m / r',labelpad=11)
    axes[0].set_ylabel('Persistence per step: (1 + m) / 2',labelpad=12)
    fig.suptitle('Persistence versus heritability: two interpretations of the intercept',
                 x=.075,y=.965,ha='left',fontsize=17,fontweight='bold',color='#213244')
    fig.text(.075,.91,'Book Table 4.5 · Eight phenotypes · Unweighted OLS of the log correlations',
             ha='left',fontsize=11,color='#536271')
    legend=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#286c91',markeredgecolor='white',markersize=8,label='Recovered parameters within bounds'),
            Line2D([0],[0],marker='X',color='none',markerfacecolor='#b7463e',markeredgecolor='white',markersize=9,label='At least one recovered parameter exceeds 1')]
    fig.legend(handles=legend,loc='lower left',bbox_to_anchor=(.07,.125),ncol=2,frameon=False,fontsize=10)
    fig.text(.075,.085,'Same OLS fits in both panels; only the horizontal-axis interpretation changes. Shading marks h² > 1.',fontsize=9,color='#536271')
    fig.text(.075,.05,'Company director: recovered r = 1.437 despite h² < 1. Panel A represents a maintained hypothesis about reporting, not verified author code.',fontsize=9,color='#536271')
    stem=OUT/'book_table45_heritability_comparison'
    for ext in ['png','pdf','svg']:fig.savefig(stem.with_suffix('.'+ext),dpi=180,facecolor='white')
    plt.close(fig)
    with (HERE/'book_table45_heritability_comparison.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['phenotype','birth_period','persistence','exp_a_theta_h2','recovered_h2','r','theta','bound_violations'])
        for f in fits:writer.writerow([f['outcome'],f['cohort'],f['b'],f['theta_h2'],f['h2'],f['r'],f['theta'],';'.join(f['bound_violations'])])
    print(stem.with_suffix('.png'))

if __name__=='__main__':main()
