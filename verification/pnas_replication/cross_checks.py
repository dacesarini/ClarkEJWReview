"""Cross-outcome checks, rounding sensitivity, and leave-one-kinship-out diagnostics."""
import json
from pathlib import Path
import sys
import numpy as np
from estimate import fit,published
from build_inputs import HERE,PUBLISHED,DISTANCE,LINEAL
sys.path.insert(0,str(HERE.parent))
from fit_all_phenotypes import CORRELATIONS,OUTCOMES

# PNAS SI Table S1, columns reordered to main Table 2 (literacy last).
COUNTS=np.array([
 [7840,7840,7895,17856,10585,15063,10772,15806,5856],
 [10678,10678,10678,17774,8530,15695,8337,15394,7467],
 [9101,9101,9160,20573,12300,35421,12030,34874,7081],
 [2285,2285,2314,4801,4734,16031,4619,15661,4450],
 [8523,8523,8523,18475,14941,27022,14635,26829,4831],
 [10843,10904,10904,25887,12744,52557,12319,51001,6215],
 [10083,10087,10165,21225,16531,35135,15935,33813,3289],
 [15708,15708,15856,33422,12843,59134,12212,57201,3645],
 [15570,15641,15761,26657,16493,35626,15349,33900,1778],
 [17843,17843,17911,37889,13533,57299,12503,53447,2830],
 [9015,9015,9015,13374,17382,33049,15953,30104,1600]])
TARGETS=np.array([[.769,.538,.466,.973],[.798,.596,.412,.978],[.732,.463,.358,.989],
 [.783,.567,.191,.718],[.812,.624,.722,.964],[.786,.572,.652,.931],
 [.774,.548,.627,.983],[.766,.531,.420,.886],[.847,.694,.427,.625]])

def main():
    checks=[]
    for i,(outcome,cohort) in enumerate(OUTCOMES):
        data=dict(rho=CORRELATIONS[:,i],N=COUNTS[:,i],n=DISTANCE,d_lin=LINEAL)
        for method in ['equal','N','normal_corr','normal_log','student_corr','student_log','inverse_se_log']:
            for iterate in ([False,True] if method.startswith('normal') else [False]):
                result=fit(data,method,iterate=iterate)
                checks.append(dict(outcome=outcome,cohort=cohort,target=dict(zip(['b','m','h2','R2'],TARGETS[i])),**result))
    rng=np.random.default_rng(20260915);rounding=[];loo=[]
    for c in PUBLISHED:
        original=published(c)
        for method in ['equal','normal_corr','normal_log']:
            vals=[]
            for _ in range(2000):
                data={**original,'rho':np.array(original['rho'])+rng.uniform(-.0005,.0005,11)}
                r=fit(data,method);vals.append([r[k] for k in ['b','m','h2','r']])
            vals=np.array(vals)
            rounding.append(dict(cohort=c,method=method,draws=2000,seed=20260915,
                note='Random exploration of rounding intervals, NOT statistical confidence intervals or rigorous extrema',
                ranges={k:[float(vals[:,i].min()),float(vals[:,i].max())] for i,k in enumerate(['b','m','h2','r'])}))
            for omit in range(11):
                data={k:[v for i,v in enumerate(arr) if i!=omit] for k,arr in original.items()}
                loo.append(dict(cohort=c,omitted_row=omit,**fit(data,method)))
    result=dict(cross_outcome=checks,rounding_sensitivity=rounding,leave_one_out=loo)
    (HERE/'cross_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    for r in rounding:print(r['cohort'],r['method'],r['ranges'])
    print('Saved cross_checks.json')

if __name__=='__main__':main()
