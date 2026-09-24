"""Theta=1 comparison: unrestricted A3 versus direct (r,m) constrained fit."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE.parents[1]/'figures'

def main():
    data=json.loads((HERE/'book_table45_ols_results.json').read_text(encoding='utf-8'))
    previous=json.loads((HERE/'book_table45_theta1_constrained_results.json').read_text(encoding='utf-8'))['fits']
    n=np.array(data['n']);d=np.array(data['d_lin']);results=[]
    for index,row in enumerate(data['fits']):
        y=np.log(row['observed'])
        def residual(z):
            r,m=np.exp(z)
            return np.log(m/r)+(n-d)*np.log((1+m)/2)+d*np.log((1+r)/2)-y
        fits=[least_squares(residual,np.log([r,m]),ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=5000)
              for r in [.2,.8,2,5] for m in [.2,.7,1.2]]
        best=min(fits,key=lambda fit:fit.fun@fit.fun)
        global_fit=differential_evolution(lambda z:np.sum(residual(z)**2),[(-10,4),(-10,4)],seed=2026,tol=1e-10)
        assert best.success
        np.testing.assert_allclose(best.fun@best.fun,global_fit.fun,atol=1e-9)
        r,m=np.exp(best.x);h2=m/r
        np.testing.assert_allclose([r,m,h2],[previous[index][k] for k in ['r','m','h2']],rtol=1e-5)
        results.append(dict(outcome=row['outcome'],cohort=row['cohort'],
                            unrestricted=dict(h2=row['theta_h2'],m=row['m'],r=row['r'],persistence=row['b']),
                            constrained=dict(h2=float(h2),m=float(m),r=float(r),persistence=float((1+m)/2))))
    (HERE/'theta1_unrestricted_vs_rm.json').write_text(json.dumps(dict(
        method='Both theta=1, equal-weight log residuals. Constrained fit directly estimates positive r and m, with h2=m/r; no upper bounds.',fits=results),indent=2)+'\n',encoding='utf-8')
    labels=['House value','IMD','Director','Occupation 1780','Occupation 1860','Education 1780','Education 1860','Literacy']
    offsets=[(-8,24),(-8,-17),(0,20),(-7,17),(8,17),(0,-22),(-12,-25),(8,17)]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(13,6),sharex=True,sharey=True)
    fig.subplots_adjust(left=.075,right=.985,top=.78,bottom=.23,wspace=.15)
    for ax,key,title in zip(axes,['unrestricted','constrained'],['1. Constraint not imposed: h² = exp(a)','2. Constraint imposed: h² = m/r']):
        for row,label,(dx,dy) in zip(results,labels,offsets):
            v=row[key];ax.scatter(v['h2'],v['persistence'],s=55,color='#27678b',edgecolors='white',linewidths=.6,zorder=3)
            ax.annotate(label,(v['h2'],v['persistence']),xytext=(dx,dy),textcoords='offset points',
                        ha='right' if dx<0 else ('center' if dx==0 else 'left'),va='center',fontsize=9,
                        arrowprops=dict(arrowstyle='-',color='#a1aab0',lw=.6))
        ax.set_xlim(0,1.04);ax.set_ylim(0,1)
        ax.set_xticks([0,.2,.4,.6,.8,1]);ax.set_yticks([0,.2,.4,.6,.8,1])
        ax.grid(color='#e2e6e9',lw=.7);ax.spines[['top','right']].set_visible(False)
        ax.set_title(title,loc='left',fontsize=12,pad=15,fontweight='bold')
        ax.set_xlabel('Heritability h²',labelpad=10)
    axes[0].set_ylabel('Persistence: (1 + m) / 2',labelpad=10)
    fig.suptitle('Both models assume θ = 1',x=.075,y=.96,ha='left',fontsize=19,fontweight='bold')
    fig.text(.075,.885,'Book Table 4.5 · Equal-weight least squares in log correlations · Same axes',fontsize=12,color='#536271')
    fig.text(.075,.09,'Left: three free coefficients. Right: directly fitted in r and m, with h² = m/r.',fontsize=10,color='#536271')
    fig.text(.075,.045,'Only the equality constraint is imposed; no upper bounds. The right-hand fit implies r > 1 for six outcomes.',fontsize=10,color='#536271')
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'book_theta1_unrestricted_vs_rm.{ext}',dpi=170,facecolor='white')
    plt.close(fig)
    print('Direct (r,m) fits verified against previous constrained fits; saved comparison.')

if __name__=='__main__':main()
