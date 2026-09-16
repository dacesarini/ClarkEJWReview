"""Fits to rounded PNAS Table 2 education correlations; no formal inference.

Requires NumPy/SciPy. SI Table S1 provides pair-count sensitivity weights.
No spouse observations are fitted. Run from any working directory.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution, least_squares

RELATIVES = ['Full siblings', 'Parent-child', 'Aunt/uncle', 'Grandparent',
             'First cousins', 'First cousins removed', 'Second cousins',
             'Second cousins removed', 'Third cousins', 'Third cousins removed',
             'Fourth cousins']
DATA = {
    '1780-1859': {
        'correlations': [.479,.538,.382,.381,.291,.234,.176,.032,.071,.094,.080],
        'pairs': [10772,8337,12030,4619,14635,12319,15935,12212,15349,12503,15953]},
    '1860-1919': {
        'correlations': [.326,.374,.228,.249,.173,.162,.186,.125,.085,.032,.014],
        'pairs': [15806,15394,34874,15661,26829,51001,33813,57201,33900,53447,30104]}}


def moments(theta, h2, m, model):
    k = (1+m)/2
    result = theta*h2*k**np.array([1,1,2,2,3,4,5,6,7,8,9])
    if model == 'phenotypic':
        result[1] = theta*(h2+m)/2
        result[3] = theta*(h2+m)*k/2
    return result


def fit(cohort, model, free_theta=False, weighted=False):
    y = np.array(DATA[cohort]['correlations'])
    w = np.array(DATA[cohort]['pairs'], dtype=float) if weighted else np.ones(11)
    w /= w.mean()

    def unpack(x):
        theta, h2, last = x if free_theta else (1.0, *x)
        # For phenotypic assortment last is r: enforce m=h2*r, 0<=r<=1.
        return theta, h2, h2*last if model == 'phenotypic' else last

    def residual(x):
        return np.sqrt(w)*(moments(*unpack(x), model)-y)

    dim = 3 if free_theta else 2
    starts = np.array(np.meshgrid(*[(.1,.5,.9)]*dim)).T.reshape(-1,dim)
    fits = [least_squares(residual, x, bounds=(0,1), ftol=1e-14,
                         xtol=1e-14, gtol=1e-14, max_nfev=5000) for x in starts]
    best = min(fits, key=lambda f: np.dot(f.fun,f.fun))
    check = differential_evolution(lambda x: np.dot(residual(x),residual(x)),
                                   [(0,1)]*dim, seed=2026, tol=1e-12, polish=True)
    if not best.success or abs(check.fun-np.dot(best.fun,best.fun)) > 1e-10:
        raise RuntimeError('Optimization verification failed')
    theta,h2,m = map(float,unpack(best.x))
    r = m/h2 if model == 'phenotypic' else h2*m
    prediction = moments(theta,h2,m,model)
    return dict(cohort=cohort, model=model, free_theta=free_theta,
                weighting='pair counts' if weighted else 'unweighted',
                theta=theta,h2=h2,m=m,r_true=r,r_observed=theta*r,
                theta_h2=theta*h2,raw_sse=float(np.sum((prediction-y)**2)),
                weighted_objective=float(np.dot(best.fun,best.fun)),
                predictions=prediction.tolist(),observed=y.tolist())


def main():
    results = [fit(cohort,model,weighted=weighted) for cohort in DATA
               for model in ('genetic','phenotypic') for weighted in (False,True)]
    results.append(fit('1860-1919','phenotypic',free_theta=True))
    genetic = next(r for r in results if r['cohort']=='1860-1919'
                   and r['model']=='genetic' and r['weighting']=='unweighted')
    product,m = genetic['theta_h2'],genetic['m']
    for theta in (product,.6,.8,1):
        np.testing.assert_allclose(moments(theta,product/theta,m,'genetic'),
                                   genetic['predictions'],rtol=1e-12,atol=1e-12)
    report = dict(source='Clark (2023), PNAS Table 2; SI Table S1 pair counts',
                  relatives=RELATIVES,inputs=DATA,fits=results,
                  genetic_free_theta=dict(theta_h2=product,m=m,
                      theta_range=[product,1],h2='theta_h2 / theta',
                      true_r_range=[product*m,m],observed_r=product*m,
                      raw_sse=genetic['raw_sse']),
                  limitations='Rounded moments; no spouse observation; equal-sex assumptions; no formal inference.')
    destination = Path(__file__).with_name('assortment_fit_results.json')
    destination.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for result in results:
        print(result['cohort'],result['model'],result['weighting'],
              'theta free' if result['free_theta'] else 'theta fixed',
              {k:round(result[k],8) for k in ('theta','h2','m','r_true','raw_sse')})
    print('Verified genetic nonidentification; saved',destination)


if __name__ == '__main__':
    main()
