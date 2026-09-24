"""Paired spouse-correlation estimates for the two theta=1 fits."""
import json
from pathlib import Path
import numpy as np
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
    fig,ax=plt.subplots(figsize=(11.5,7.7))
    fig.subplots_adjust(left=.27,right=.95,top=.77,bottom=.19)
    ax.axvspan(1,4.5,color='#fbefed',zorder=0)
    ax.axvline(1,color='#a5433c',lw=1.3,ls=(0,(4,3)),zorder=1)
    ax.set_xlim(0,4.5);ax.set_ylim(-.6,7.8)
    for i,(f,marker,color) in enumerate(zip(fits,markers,colors)):
        y=7-i;u=f['unrestricted']['r'];c=f['constrained']['r']
        ax.hlines(y,0,4.5,color='#e6e9ec',lw=.7,zorder=0)
        ax.plot([u,c],[y,y],color=color,lw=2.5,alpha=.65,zorder=2)
        ax.scatter(u,y,marker=marker,s=125,facecolors='white',edgecolors=color,linewidths=1.8,zorder=4)
        ax.scatter(c,y,marker=marker,s=48,color='black',linewidths=.6,zorder=5)
        ax.annotate(f'{u:.3f}',(u,y),xytext=(-8,10),textcoords='offset points',ha='right',va='bottom',fontsize=10,color=color,fontweight='bold')
        ax.annotate(f'{c:.3f}',(c,y),xytext=(8,-10),textcoords='offset points',ha='left',va='top',fontsize=10,color='#222222',fontweight='bold')
    ax.set_yticks(range(8),labels=labels[::-1]);ax.tick_params(axis='y',length=0,pad=14)
    ax.set_xticks(np.arange(0,4.51,.5));ax.tick_params(axis='x',length=4,color='#6c7780')
    ax.set_xlabel('Spouse phenotypic correlation r',labelpad=13)
    ax.spines[['top','right','left']].set_visible(False)
    ax.text(1.04,7.57,'Above the maximum possible correlation',color='#a5433c',fontsize=10,va='center')
    fig.suptitle('Similar heritabilities. Very different spouse correlations.',x=.06,y=.96,ha='left',fontsize=18,fontweight='bold')
    fig.text(.06,.898,'Both θ = 1 · Book Table 4.5 · Equal-weight least squares in log correlations',fontsize=11,color='#536271')
    fig.text(.06,.842,'Estimates above 1:  unconstrained 1 of 8  →  constrained 6 of 8',fontsize=13,fontweight='bold',color='#a5433c')
    handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor='white',markeredgecolor='#0072B2',markeredgewidth=1.7,markersize=9,label='Unconstrained (colored outline)'),
             Line2D([0],[0],marker='o',color='none',markerfacecolor='black',markeredgecolor='black',markersize=6,label='Constrained (solid black)')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.055,.072),ncol=2,frameon=False,fontsize=10)
    fig.text(.06,.045,'Lines connect the same trait. Constrained fit imposes m = h²r; no upper bounds are imposed.',fontsize=10,color='#536271')
    fig.text(.06,.017,'Point estimates only. Colors and shapes match the heritability–persistence comparison.',fontsize=9,color='#536271')
    OUT.mkdir(exist_ok=True)
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'book_theta1_spouse_r_comparison.{ext}',dpi=180,facecolor='white')
    plt.close(fig)

if __name__=='__main__':main()
