"""Diagnostic fits allowing h2>1, retaining PH and valid m,r,theta.

Parameterize by q=theta*h2, m, r, then h2=m/r and theta=q*r/m.
Linear constraints r<=m/q enforce theta<=1 for fixed q,m; use SLSQP
with constraints q*r<=m and bounded m,r. Bounds on q from sibling
prediction and baseline SSE are derived per outcome, not arbitrary.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize
from fit_all_phenotypes import OUTCOMES, CORRELATIONS


def predicted(x):
    q,m,r=x
    k=(1+m)/2
    p=q*k**np.array([1,1,2,2,3,4,5,6,7,8,9])
    p[1]=q*(1+r)/2
    p[3]=q*k*(1+r)/2
    return p


def main():
    root=Path(__file__).resolve().parent
    original=json.loads((root/'all_phenotypes_free_theta_results.json').read_text(encoding='utf-8'))
    results=[]
    for column,(name,cohort) in enumerate(OUTCOMES):
        y=CORRELATIONS[:,column]
        baseline=original['fits'][column]
        # Any improvement must satisfy |q*(1+m)/2-y_sib| <= sqrt(SSE).
        qmax=2*(y[0]+np.sqrt(baseline['raw_sse']))+1e-6
        def residual(x):return predicted(x)-y
        def objective(x):return np.sum(residual(x)**2)
        bounds=([1e-10,1e-10,0],[qmax,1,1])
        starts=[[q,m,r] for q in (.15,min(.6,qmax*.9)) for m in (.2,.6,.9) for r in (0,.4,.9)]
        unconstrained=[least_squares(residual,x,bounds=bounds,ftol=1e-14,xtol=1e-14,gtol=1e-14,max_nfev=5000) for x in starts]
        best=min(unconstrained,key=lambda f:np.dot(f.fun,f.fun))
        q,m,r=best.x
        # If the box optimum satisfies theta<=1, it is also the optimum
        # with that extra restriction. Otherwise optimize it explicitly.
        if q*r > m+1e-9:
            starts.append([baseline['theta_h2'],baseline['m'],baseline['r_true']])
            fits=[minimize(objective,x,method='SLSQP',bounds=list(zip(*bounds)),
                           constraints=[{'type':'ineq','fun':lambda x:x[1]-x[0]*x[2]}],
                           options={'ftol':1e-13,'maxiter':3000}) for x in starts]
            feasible=[f for f in fits if f.success and f.x[1]-f.x[0]*f.x[2]>=-1e-8]
            best=min(feasible,key=lambda f:f.fun)
            q,m,r=best.x
        else:
            # Independently check with a different optimizer.
            check=minimize(objective,best.x,method='SLSQP',bounds=list(zip(*bounds)),
                           constraints=[{'type':'ineq','fun':lambda x:x[1]-x[0]*x[2]}],
                           options={'ftol':1e-13,'maxiter':3000})
            assert abs(check.fun-objective(best.x))<1e-9
        boundary=r<1e-7
        h2=None if boundary else float(m/r)
        theta=0.0 if boundary else float(q*r/m)
        result=dict(outcome=name,cohort=cohort,theta_h2=float(q),m=float(m),r_true=float(r),
                    h2=h2,theta=theta,infinite_h2_limit=bool(boundary),
                    r_observed=float(q*r*r/m),raw_sse=float(objective(best.x)),
                    bounded_sse=baseline['raw_sse'],predictions=predicted(best.x).tolist())
        assert result['raw_sse']<=baseline['raw_sse']+1e-9
        assert q*r<=m+1e-8
        results.append(result)
        print(name,cohort, {k:result[k] for k in ('h2','m','theta','r_true','theta_h2','raw_sse','infinite_h2_limit')})
    out=root/'relaxed_heritability_results.json'
    out.write_text(json.dumps(dict(method='Unweighted correlation-level PH fit; h2 upper bound removed, theta<=1 and 0<=m,r<=1 retained. q=theta*h2. Zero-r optima are limiting fits with h2 tending to infinity.',fits=results),indent=2)+'\n',encoding='utf-8')
    print('Saved',out)


if __name__=='__main__':
    main()
