"""Unrestricted A3 OLS for book Table 4.5, then free-theta PH recovery.

No admissibility constraints are imposed on OLS. Invalid structural recoveries
are reported, never clipped. Uses book inputs, NOT final PNAS Table 2.
"""
import json
import re
import hashlib
import numpy as np
from build_inputs import ROOT,HERE,DISTANCE,LINEAL

def main():
    source=ROOT/'references/text/Clark_Book_2026.txt'
    text=source.read_text(encoding='utf-8')
    table=text[text.index('Table 4.5:'):];table=table[:table.index('Unrelated')]
    lines=[l for l in table.splitlines() if re.match(r'^\s*(Full Sibling|Child\s|Sibling-rem|Grandchild|Cousin)',l)]
    labels=[re.split(r'\s+0\.',l.strip(),maxsplit=1)[0] for l in lines]
    y=np.array([[float(v) for v in re.findall(r'0\.\d+',l)] for l in lines])
    assert y.shape==(11,8),y.shape
    outcomes=[('House value','1910-1996'),('IMD','1910-1996'),('Company director','1910-1996'),
              ('Occupational status','1780-1859'),('Occupational status','1860-1919'),
              ('Higher education','1780-1859'),('Higher education','1860-1919'),('Literacy','1725-1869')]
    X=np.column_stack([np.ones(11),DISTANCE,LINEAL]);coef=np.linalg.lstsq(X,np.log(y),rcond=None)[0]
    results=[]
    for j,(outcome,cohort) in enumerate(outcomes):
        a,beta,c=map(float,coef[:,j]);b=float(np.exp(beta));m=2*b-1;r=2*b*np.exp(c)-1
        product=float(np.exp(a));h2=m/r;theta=product/h2
        fitted=X@coef[:,j];errors=np.log(y[:,j])-fitted
        violations=[]
        for name,value in [('h2',h2),('theta',theta),('m',m),('r',r)]:
            if not 0<=value<=1:violations.append(name)
        reconstructed=theta*h2*((1+m)/2)**np.array(DISTANCE)*((1+r)/(1+m))**np.array(LINEAL)
        np.testing.assert_allclose(reconstructed,np.exp(fitted),rtol=1e-12)
        np.testing.assert_allclose(X.T@errors,np.zeros(3),atol=1e-12)
        result=dict(outcome=outcome,cohort=cohort,a=a,beta=beta,c=c,b=b,m=m,r=float(r),
                    theta_h2=product,h2=float(h2),theta=float(theta),bound_violations=violations,
                    R2=float(1-errors@errors/np.sum((np.log(y[:,j])-np.log(y[:,j]).mean())**2)),
                    observed=y[:,j].tolist(),predicted=np.exp(fitted).tolist())
        results.append(result)
        print(outcome,cohort,{k:round(result[k],6) for k in ['m','r','theta_h2','h2','theta']},violations)
    report=dict(source='Book Table 4.5, printed p.54; A1-A3, pp.296-297',
                extract_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),relatives=labels,
                n=DISTANCE,d_lin=LINEAL,method='Equal-weight OLS on log correlations; unrelated excluded; theta free PH transformation, no bounds enforced',
                transformations=dict(b='exp(beta)',m='2*b-1',r='(1+m)*exp(c)-1',theta_h2='exp(a)',h2='m/r',theta='exp(a)/h2'),fits=results)
    (HERE/'book_table45_ols_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
