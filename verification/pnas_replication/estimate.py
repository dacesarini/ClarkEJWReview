"""Reproducible Clark Eq.3 sensitivity engine; candidate methods, not author code.

Usage: python estimate.py --suite
       python estimate.py --cohort 1860-1919 --weight normal_log --model ph_fixed
       python estimate.py --csv custom.csv --weight supplied --model unrestricted
CSV columns: rho,N,n,d_lin; optional se_binary,weight. Rows can be omitted via --omit.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from build_inputs import PUBLISHED, DISTANCE, LINEAL, LABELS

HERE=Path(__file__).resolve().parent
WEIGHTS=['equal','N','normal_corr','normal_log','student_corr','student_log',
         'binary_corr','binary_log','inverse_se_corr','inverse_se_log','N_rho2','supplied']
MODELS=['unrestricted','ph_fixed','ph_free','genetic']

def weights(kind,rho,N,se=None,supplied=None):
    if kind=='equal':return np.ones(len(rho))
    if kind=='N':return N.copy()
    if kind=='supplied':return np.asarray(supplied,dtype=float)
    if kind=='N_rho2':return N*rho**2
    if kind.startswith('binary'):
        if se is None or not np.all(np.isfinite(se)):raise ValueError('Binary SEs unavailable; refusing silent substitution')
        var=np.asarray(se)**2
    elif kind.startswith('student'):var=(1-rho**2)/(N-2)
    else:var=(1-rho**2)**2/(N-3)
    if kind.endswith('log'):var=var/rho**2
    return 1/np.sqrt(var) if kind.startswith('inverse_se') else 1/var

def log_prediction(p,n,d,model):
    if model=='unrestricted':
        a,beta,c=p
    elif model in ('ph_fixed','ph_free'):
        theta,h2,r=p if model=='ph_free' else (1.,*p)
        m=h2*r;a=np.log(theta*h2);beta=np.log((1+m)/2);c=np.log((1+r)/(1+m))
    else:
        h2,m=p;a=np.log(h2);beta=np.log((1+m)/2);c=0.
    return a+n*beta+d*c

def solve(rho,n,d,w,model,scale):
    X=np.column_stack([np.ones(len(n)),n,d]);root=np.sqrt(w/w.mean())
    if model=='unrestricted' and scale=='log':
        if np.linalg.matrix_rank(X)<3:raise ValueError('Regressors do not identify all three coefficients')
        p=np.linalg.lstsq(X*root[:,None],np.log(rho)*root,rcond=None)[0]
        return p,True
    def residual(p):
        pred=log_prediction(p,n,d,model)
        return root*(pred-np.log(rho) if scale=='log' else np.exp(pred)-rho)
    if model=='unrestricted':
        starts=[np.linalg.lstsq(X,np.log(rho),rcond=None)[0]]
        # Broad numerical bounds, not structural correlation restrictions.
        bounds=([-15,-5,-5],[3,1,5])
    else:
        dim=3 if model=='ph_free' else 2
        starts=[np.full(dim,.3),np.full(dim,.7),np.array([.9,.5,.8]) if dim==3 else np.array([.5,.9])]
        bounds=(np.full(dim,1e-8),np.ones(dim))
    results=[least_squares(residual,start,bounds=bounds,ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=5000) for start in starts]
    best=min(results,key=lambda f:float(f.fun@f.fun))
    if not best.success:raise RuntimeError(best.message)
    return best.x,best.success

def fit(data,weight='equal',model='unrestricted',scale='log',iterate=False):
    rho=np.asarray(data['rho'],float);N=np.asarray(data['N'],float)
    n=np.asarray(data['n'],float);d=np.asarray(data['d_lin'],float)
    if not (rho.shape==N.shape==n.shape==d.shape):raise ValueError('Input vectors have different lengths')
    design=np.column_stack([np.ones(len(n)),n,d]) if model!='genetic' else np.column_stack([np.ones(len(n)),n])
    if model in ('unrestricted','ph_free','genetic') and np.linalg.matrix_rank(design)<design.shape[1]:
        raise ValueError('Selected moments do not identify this model')
    if np.any((rho<=0)|(rho>=1)):raise ValueError('Positive correlations below one required for this Eq.3 implementation')
    if np.any(N<=3):raise ValueError('N must exceed three')
    if iterate and weight not in ('normal_corr','normal_log','student_corr','student_log','N_rho2'):
        raise ValueError('Iteration supported only for explicitly model-based normal/student/N*rho^2 weights')
    current=rho.copy();previous=None;converged=not iterate
    for iteration in range(1,201 if iterate else 2):
        w=weights(weight,current,N,data.get('se_binary'),data.get('weight'))
        if w.shape!=rho.shape or not np.all(np.isfinite(w)) or np.any(w<=0):raise ValueError('Invalid weights')
        p,success=solve(rho,n,d,w,model,scale);logpred=log_prediction(p,n,d,model);pred=np.exp(logpred)
        if not iterate:break
        if previous is not None and np.max(np.abs(pred-previous))<1e-11:
            converged=True;break
        if np.any((pred<=0)|(pred>=1)):raise ValueError('IRLS predicted correlations outside (0,1)')
        previous=pred.copy();current=pred
    if model=='unrestricted':
        a,beta,c=map(float,p);h2=float(np.exp(a));m=float(2*np.exp(beta)-1);r=float((1+m)*np.exp(c)-1);theta=1.
    elif model in ('ph_fixed','ph_free'):
        theta,h2,r=map(float,p) if model=='ph_free' else (1.,*map(float,p));m=h2*r
        a=float(np.log(theta*h2));beta=float(np.log((1+m)/2));c=float(np.log((1+r)/(1+m)))
    else:
        h2,m=map(float,p);theta=1.;r=h2*m;a=float(np.log(h2));beta=float(np.log((1+m)/2));c=0.
    response=np.log(rho) if scale=='log' else rho;fitted=logpred if scale=='log' else pred
    R2=1-np.sum(w*(response-fitted)**2)/np.sum(w*(response-np.average(response,weights=w))**2)
    return dict(weight=weight,model=model,scale=scale,iterate=iterate,converged=converged,iterations=iteration,
                a=a,beta=beta,c=c,b=float(np.exp(beta)),h2=h2,m=m,r=r,theta=theta,
                exp_a=float(np.exp(a)),constraint_gap_ph=m-h2*r,weighted_R2=float(R2),
                weights_normalized=(w/w.mean()).tolist(),observed=rho.tolist(),predictions=pred.tolist(),
                objective_normalized=float(np.sum((w/w.mean())*(response-fitted)**2)),
                structural_bounds_ok=bool(0<h2<=1 and 0<=m<=1 and 0<=r<=1 and 0<theta<=1))

def published(cohort):
    return dict(rho=PUBLISHED[cohort]['rho'],N=PUBLISHED[cohort]['N'],n=DISTANCE,d_lin=LINEAL)

def suite():
    audit=json.loads((HERE/'inputs_audit.json').read_text(encoding='utf-8'))
    output=[];failures=[]
    for cohort in PUBLISHED:
        datasets={'published':published(cohort)}
        for rule in ('flags','birth_rule'):
            rows={r['raw_label']:r for r in audit['pair_audit']['group_stats'] if r['cohort']==cohort and r['rule']==rule}
            datasets['literal_raw_'+rule]=dict(rho=[rows[l]['rho'] for l in LABELS],N=[rows[l]['N'] for l in LABELS],
                se_binary=[rows[l]['se_binary'] for l in LABELS],n=DISTANCE,d_lin=LINEAL)
        for source,data in datasets.items():
            for weight in WEIGHTS:
                if weight=='supplied' or (weight.startswith('binary') and source=='published'):continue
                for model in MODELS:
                    for scale in ('log','level'):
                        recipes=[False,True] if weight in ('normal_corr','normal_log') else [False]
                        for iterate in recipes:
                            tag=dict(cohort=cohort,source=source,weight=weight,model=model,scale=scale,iterate=iterate)
                            try:
                                result=fit(data,weight,model,scale,iterate)
                                result.update(cohort=cohort,source=source)
                                target=PUBLISHED[cohort]['clark']
                                result['differences_from_clark']={k:(result['exp_a'] if k=='h2' else result[k])-target[k] for k in ('b','m','h2')}
                                result['comparison_note']='Clark reported h2 is compared with exp(a), not latent h2 in the free-theta model.'
                                output.append(result)
                            except (ValueError,RuntimeError) as error:failures.append(dict(**tag,error=str(error)))
    report=dict(note='Candidate procedures only; no exact author-code replication. Literal raw labels are known not to match all published inputs. All diagonal variance formulas ignore dependence between relatives.',fits=output,failures=failures)
    (HERE/'fit_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    keys=['cohort','source','model','scale','weight','iterate','converged','a','beta','c','b','m','h2','r','theta','weighted_R2']
    with (HERE/'fit_summary.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');writer.writeheader();writer.writerows(output)
    for cohort in PUBLISHED:
        print('\nClosest candidate fits to printed b and h2:',cohort)
        candidates=[f for f in output if f['cohort']==cohort and f['source']=='published' and f['model']=='unrestricted' and f['converged']]
        for f in sorted(candidates,key=lambda f:sum(f['differences_from_clark'][k]**2 for k in ('b','h2')))[:6]:
            print(f['scale'],f['weight'],'iterated',f['iterate'],{k:round(f[k],6) for k in ('b','m','h2','r','weighted_R2')})
    print(f'Saved {len(output)} fits; {len(failures)} failed candidates explicitly recorded.')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite',action='store_true');parser.add_argument('--cohort',choices=list(PUBLISHED),default='1780-1859')
    parser.add_argument('--csv',type=Path);parser.add_argument('--weight',choices=WEIGHTS,default='equal')
    parser.add_argument('--model',choices=MODELS,default='unrestricted');parser.add_argument('--scale',choices=['log','level'],default='log')
    parser.add_argument('--iterate',action='store_true');parser.add_argument('--omit',type=int,nargs='*',default=[],help='Zero-based row indices')
    args=parser.parse_args()
    if args.suite:suite();return
    if args.csv:
        with args.csv.open(encoding='utf-8-sig') as f:rs=list(csv.DictReader(f))
        data={key:[float(r[key]) for r in rs] for key in rs[0] if key in ('rho','N','n','d_lin','se_binary','weight')}
    else:data=published(args.cohort)
    data={k:[v for i,v in enumerate(values) if i not in args.omit] for k,values in data.items()}
    print(json.dumps(fit(data,args.weight,args.model,args.scale,args.iterate),indent=2))

if __name__=='__main__':main()
