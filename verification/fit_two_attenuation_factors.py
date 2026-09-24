"""Phenotypic assortment with separate same/cross-generation attenuation.

Pair-specific attenuation extension, not automatically an independent
individual measurement-error construction. No spouse moments fitted.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, differential_evolution
from fit_all_phenotypes import OUTCOMES, CORRELATIONS
from fit_assortment_models import RELATIVES, moments

SAME_GENERATION = np.array([True,False,False,False,True,False,True,False,True,False,True])


def prediction(x):
    h2,r,theta_same,theta_cross=x
    m=h2*r
    return moments(1,h2,m,'phenotypic')*np.where(SAME_GENERATION,theta_same,theta_cross)


def main():
    root=Path(__file__).resolve().parent
    old=json.loads((root/'all_phenotypes_free_theta_results.json').read_text(encoding='utf-8'))
    results=[]
    for column,(outcome,cohort) in enumerate(OUTCOMES):
        y=CORRELATIONS[:,column]
        def residual(x):return prediction(x)-y
        starts=np.array(np.meshgrid(*[(.2,.7)]*4)).T.reshape(-1,4)
        fits=[least_squares(residual,x,bounds=(0,1),ftol=1e-14,xtol=1e-14,
                            gtol=1e-14,max_nfev=5000) for x in starts]
        best=min(fits,key=lambda f:np.dot(f.fun,f.fun))
        global_fit=differential_evolution(lambda x:np.sum(residual(x)**2),[(0,1)]*4,
                                          seed=2026,tol=1e-12,polish=True)
        assert best.success and abs(global_fit.fun-np.sum(best.fun**2))<1e-9
        h2,r,ts,tc=map(float,best.x)
        result=dict(outcome=outcome,cohort=cohort,h2=h2,m=h2*r,r_true=r,
                    theta_same=ts,theta_cross=tc,raw_sse=float(np.sum(best.fun**2)),
                    single_theta_sse=old['fits'][column]['raw_sse'],
                    observed=y.tolist(),predictions=prediction(best.x).tolist())
        result['boundaries']=[key for key in ('h2','r_true','theta_same','theta_cross')
                              if min(result[key],1-result[key])<1e-7]
        assert result['raw_sse']<=result['single_theta_sse']+1e-9
        results.append(result)
        print(outcome,cohort,{key:round(result[key],6) for key in
             ('h2','m','theta_same','theta_cross','r_true','raw_sse')},result['boundaries'])
    report=dict(source='PNAS Table 2 rounded correlations',
                method='Phenotypic assortment, m=h2*r; h2,r,theta_same,theta_cross in [0,1]; unweighted correlation-level least squares; multistart and differential-evolution verification.',
                relatives=RELATIVES,same_generation=SAME_GENERATION.tolist(),fits=results,
                qualification='Two pair-class attenuation factors require additional measurement assumptions; they are not automatically compatible with a single classical independent-error model for individuals. No standard errors or spouse moments included.')
    path=root/'two_attenuation_factors_results.json'
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Saved',path)


if __name__=='__main__':
    main()
