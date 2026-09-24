"""Paired overlay of unrestricted and equality-constrained theta=1 estimates."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent
OUT=HERE.parents[1]/'figures'

def main():
    fits=json.loads((HERE/'theta1_unrestricted_vs_rm.json').read_text(encoding='utf-8'))['fits']
    labels=['House value','IMD','Company director','Occupation 1780–1859',
            'Occupation 1860–1919','Education 1780–1859','Education 1860–1919','Literacy']
    markers=['o','s','^','D','v','P','X','*']
    colors=['#0072B2','#009E73','#D55E00','#CC79A7','#56B4E9','#E69F00','#6A51A3','#8C564B']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(11,6.7))
    fig.subplots_adjust(left=.10,right=.71,bottom=.18,top=.82)
    for f,marker,color in zip(fits,markers,colors):
        u,c=f['unrestricted'],f['constrained']
        ax.plot([u['h2'],c['h2']],[u['persistence'],c['persistence']],color=color,lw=1.7,alpha=.75,zorder=2)
        # Larger hollow colored symbol remains visible when estimates nearly coincide.
        ax.scatter(u['h2'],u['persistence'],marker=marker,s=125,facecolors='white',edgecolors=color,linewidths=1.8,zorder=4)
        ax.scatter(c['h2'],c['persistence'],marker=marker,s=38,color='black',linewidths=.6,zorder=5)
    ax.set_xlim(0,1.04);ax.set_ylim(0,1)
    ax.set_xticks([0,.2,.4,.6,.8,1]);ax.set_yticks([0,.2,.4,.6,.8,1])
    ax.set_xlabel('Heritability h²',labelpad=12);ax.set_ylabel('Persistence: (1 + m) / 2',labelpad=12)
    ax.grid(color='#e2e6e9',lw=.7);ax.spines[['top','right']].set_visible(False)
    trait_handles=[Line2D([0],[0],marker=m,color='none',markerfacecolor='white',markeredgecolor=c,
                           markeredgewidth=1.5,markersize=9,label=l) for m,c,l in zip(markers,colors,labels)]
    trait_legend=ax.legend(handles=trait_handles,title='Trait',loc='upper left',bbox_to_anchor=(1.025,1.015),
                           frameon=False,fontsize=10,labelspacing=1.1,handletextpad=.7)
    ax.add_artist(trait_legend)
    ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='white',markeredgecolor='#0072B2',markeredgewidth=1.7,markersize=9,label='Unconstrained (colored outline)'),
                       Line2D([0],[0],marker='o',color='none',markerfacecolor='black',markeredgecolor='black',markersize=5,label='Constrained (solid black)')],
              title='Estimation',loc='upper left',bbox_to_anchor=(1.025,.26),frameon=False,fontsize=9,labelspacing=1.0)
    fig.suptitle('How much does the constraint move the estimates?',x=.10,y=.96,ha='left',fontsize=18,fontweight='bold')
    fig.text(.10,.895,'Both θ = 1 · Book Table 4.5 · Equal-weight least squares in log correlations',fontsize=11,color='#536271')
    fig.text(.10,.07,'Lines connect estimates for the same trait. Constrained fit imposes m = h²r, with h² = m/r.',fontsize=10,color='#536271')
    fig.text(.10,.035,'No upper parameter bounds imposed; six constrained r estimates exceed 1. Point estimates only.',fontsize=10,color='#536271')
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'book_theta1_paired_overlay.{ext}',dpi=180,facecolor='white')
    plt.close(fig)

if __name__=='__main__':main()
