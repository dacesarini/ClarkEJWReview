"""Simple Figure A1 analogue: unrestricted Table 4.5 OLS, theta=1 labels."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / 'figures'

def main():
    fits=json.loads((HERE/'book_table45_ols_results.json').read_text(encoding='utf-8'))['fits']
    labels=['House value','IMD','Company director','Occupation 1780–1859',
            'Occupation 1860–1919','Education 1780–1859','Education 1860–1919','Literacy']
    offsets=[(-12,20),(-12,-18),(0,17),(-8,15),(8,14),(0,-24),(-14,-22),(8,13)]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(9,6.2))
    fig.subplots_adjust(left=.12,right=.96,bottom=.20,top=.83)
    ax.scatter([f['theta_h2'] for f in fits],[f['b'] for f in fits],s=70,
               color='#27678b',edgecolors='white',linewidths=.8,zorder=3)
    for f,label,(dx,dy) in zip(fits,labels,offsets):
        ax.annotate(label,(f['theta_h2'],f['b']),xytext=(dx,dy),textcoords='offset points',
                    ha='right' if dx<0 else ('center' if dx==0 else 'left'),va='center',fontsize=10,
                    arrowprops=dict(arrowstyle='-',color='#a1aab0',lw=.65))
    ax.set_xlim(0,1.04);ax.set_ylim(0,1)
    ax.set_xticks([0,.2,.4,.6,.8,1]);ax.set_yticks([0,.2,.4,.6,.8,1])
    ax.set_xlabel('Heritability: h² = exp(a), assuming θ = 1',labelpad=12)
    ax.set_ylabel('Persistence: (1 + m) / 2',labelpad=12)
    ax.grid(color='#e2e6e9',linewidth=.7,zorder=0)
    ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Persistence versus heritability',x=.12,y=.96,ha='left',fontsize=18,fontweight='bold')
    fig.text(.12,.895,'Our unweighted OLS estimates · Book Table 4.5',fontsize=12,color='#536271')
    fig.text(.12,.065,'Unrestricted A3: m = h²r is not imposed. Point estimates only.',fontsize=9,color='#536271')
    OUT.mkdir(exist_ok=True)
    for ext in ('png','pdf','svg'):
        fig.savefig(OUT/f'book_table45_theta1_simple.{ext}',dpi=170,facecolor='white')
    plt.close(fig)

if __name__=='__main__':main()
