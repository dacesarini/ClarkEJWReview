"""Refit Table 4.5 with theta=1 and m=h2*r, without upper bounds of one."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE.parents[1]/'figures'

def main():
    source=json.loads((HERE/'book_table45_ols_results.json').read_text(encoding='utf-8'))
    n=np.array(source['n']);d=np.array(source['d_lin']);fits=[]
    for row in source['fits']:
        y=np.log(row['observed'])
        def residual(z):
            h2,m=np.exp(z);r=m/h2
            return np.log(h2)+n*np.log((1+m)/2)+d*np.log((1+r)/(1+m))-y
        candidates=[least_squares(residual,np.log(start),ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=5000)
                    for start in [(h,m) for h in [.1,.5,1,2] for m in [.1,.5,1]]]
        best=min(candidates,key=lambda f:f.fun@f.fun)
        check=differential_evolution(lambda z:np.sum(residual(z)**2),[(-10,3),(-10,3)],seed=2026,tol=1e-10)
        assert best.success
        np.testing.assert_allclose(best.fun@best.fun,check.fun,atol=1e-9)
        h2,m=np.exp(best.x);r=m/h2
        unrestricted_sse=np.sum((np.log(row['predicted'])-y)**2)
        assert best.fun@best.fun>=unrestricted_sse-1e-10
        fitted=dict(outcome=row['outcome'],cohort=row['cohort'],h2=float(h2),m=float(m),r=float(r),
                    theta=1.,persistence=float((1+m)/2),sse_log=float(best.fun@best.fun),
                    predicted=np.exp(y+best.fun).tolist())
        fits.append(fitted)
        print(row['outcome'],row['cohort'],{k:round(fitted[k],6) for k in ['h2','m','r','persistence']})
    (HERE/'book_table45_theta1_constrained_results.json').write_text(json.dumps(dict(
        method='Equal-weight least squares in log correlations, theta=1, m=h2*r; h2 and m positive, no upper bounds of one; multistart and global optimization checked',fits=fits),indent=2)+'\n',encoding='utf-8')
    labels=['House value','IMD','Company director','Occupation 1780–1859',
            'Occupation 1860–1919','Education 1780–1859','Education 1860–1919','Literacy']
    offsets=[(-12,20),(-12,-18),(0,17),(-8,15),(8,14),(0,-24),(-14,-22),(8,13)]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(9,6.2));fig.subplots_adjust(left=.12,right=.96,bottom=.20,top=.83)
    ax.scatter([f['h2'] for f in fits],[f['persistence'] for f in fits],s=70,color='#27678b',edgecolors='white',linewidths=.8,zorder=3)
    for f,label,(dx,dy) in zip(fits,labels,offsets):
        ax.annotate(label,(f['h2'],f['persistence']),xytext=(dx,dy),textcoords='offset points',
                    ha='right' if dx<0 else ('center' if dx==0 else 'left'),va='center',fontsize=10,
                    arrowprops=dict(arrowstyle='-',color='#a1aab0',lw=.65))
    ax.set_xlim(0,max(1.04,max(f['h2'] for f in fits)+.05));ax.set_ylim(0,1)
    ax.set_xticks([0,.2,.4,.6,.8,1]);ax.set_yticks([0,.2,.4,.6,.8,1])
    ax.set_xlabel('Heritability: h² = exp(a), assuming θ = 1',labelpad=12)
    ax.set_ylabel('Persistence: (1 + m) / 2',labelpad=12)
    ax.grid(color='#e2e6e9',linewidth=.7,zorder=0);ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Persistence versus heritability',x=.12,y=.96,ha='left',fontsize=18,fontweight='bold')
    fig.text(.12,.895,'Our constrained least-squares estimates · Book Table 4.5',fontsize=12,color='#536271')
    fig.text(.12,.065,'θ = 1 and m = h²r imposed. Equal-weight log residuals; no upper parameter bounds. Point estimates only.',fontsize=9,color='#536271')
    OUT.mkdir(exist_ok=True)
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'book_table45_theta1_constrained.{ext}',dpi=170,facecolor='white')
    plt.close(fig)

if __name__=='__main__':main()
