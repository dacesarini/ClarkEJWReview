"""Checks for the review's classical-measurement-error extension."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import verify_fisher_model as base

OUT=Path(__file__).resolve().parents[1]/"_build"
OUT.mkdir(exist_ok=True)
H,r,T=sp.symbols("H r T",positive=True)
L=1-T
k=(1+H*r)/2
checks=[]
def check(name,actual,expected):
    residual=sp.simplify(actual-expected)
    assert residual==0,(name,residual)
    checks.append(name)
a,e,u,g,b,c,d=sp.symbols("a e u g b c d",real=True)
w=sp.Matrix([a,e,u])
K=sp.Matrix([[g,c,0],[d,b,0],[0,0,0]])
check("all error cross terms removed from observed covariance",
      (w.T*K*w)[0],a*a*g+a*e*(c+d)+e*e*b)
check("observed variance normalization",L*H+L*(1-H)+T,1)
check("observed additive genetic variance share",L*H,(1-T)*H)
check("latent heritability from observed genetic share",(L*H)/L,H)
hP,eP=sp.sqrt(H),sp.sqrt(1-H)
wP=sp.Matrix([sp.sqrt(L)*hP,sp.sqrt(L)*eP,sp.sqrt(T)])
Cplus=r*sp.Matrix([[H,hP*eP,0],[hP*eP,1-H,0],[0,0,0]])
check("observed spouse correlation",(wP.T*Cplus*wP)[0],L*r)
check("true genetic spouse correlation from observed quantities",
      (L*H)*(L*r)/L**2,H*r)
targets={
 "spouses":r,"parent_child":H*(1+r)/2,"full_siblings":H*k,
 "grandparent":H*(1+r)*k/2,"avuncular":H*k**2,
 "first_cousins":H*k**3,"second_cousins":H*k**5,
 "third_cousins":H*k**7,"fourth_cousins":H*k**9,
 "first_cousins_removed":H*k**4,"second_cousins_removed":H*k**6,
 "third_cousins_removed":H*k**8}
nd={
 "parent_child":(1,1),"full_siblings":(1,0),"grandparent":(2,1),
 "avuncular":(2,0),"first_cousins":(3,0),"second_cousins":(5,0),
 "third_cousins":(7,0),"fourth_cousins":(9,0),
 "first_cousins_removed":(4,0),"second_cousins_removed":(6,0),
 "third_cousins_removed":(8,0)}
alpha=sp.log(L*H)
beta=sp.log(k)
gamma=sp.log((1+r)/(1+H*r))
for name,latent in targets.items():
    observed=L*latent
    check(name+" no noise limit",observed.subs(T,0),latent)
    check(name+" all noise limit",observed.subs(T,1),0)
    if name in nd:
        n,D=nd[name]
        check(name+" log regression representation",
              sp.exp(alpha+n*beta+D*gamma),observed)
for degree in range(1,5):
    check("cousin degree "+str(degree)+" transmission recurrence",
          H*k*(k*k)**degree,H*k**(2*degree+1))
m_rec=2*sp.exp(beta)-1
r_rec=2*sp.exp(beta+gamma)-1
H_rec=sp.simplify(m_rec/r_rec)
L_rec=sp.exp(alpha)/H_rec
check("recover genetic spouse correlation",m_rec,H*r)
check("recover latent spouse correlation",r_rec,r)
check("recover latent heritability",H_rec,H)
check("recover signal share",L_rec,L)
check("recover noise share",1-L_rec,T)
J=sp.Matrix([alpha,beta,gamma]).jacobian([H,r,T])
check("identification Jacobian determinant",J.det(),
      -r/(L*(1+r)*(1+H*r)))
check("random mating loses the lineal contrast",gamma.subs(r,0),0)
check("random mating decay independent of heritability",beta.subs(r,0),sp.log(sp.Rational(1,2)))
vbar=(1+L*r)/2
pc=L*H*(1+r)/2
corrbar=pc/sp.sqrt(vbar)
slopebar=pc/vbar
check("midparent covariance as average of parent covariances",(pc+pc)/2,pc)
check("midparent correlation",corrbar,L*H*(1+r)/sp.sqrt(2*(1+L*r)))
check("midparent no noise correlation",corrbar.subs(T,0),H*sp.sqrt((1+r)/2))
check("midparent slope",slopebar,L*H*(1+r)/(1+L*r))
check("midparent no noise slope",slopebar.subs(T,0),H)
check("midparent error slope differs from lambda times h2",
      slopebar-L*H,L*H*r*T/(1+L*r))
check("observed parent minus sibling gap",
      L*(targets["parent_child"]-targets["full_siblings"]),L*H*r*(1-H)/2)
# Exact observational equivalence for all collateral moments.
p1={H:sp.Rational(9,25),r:sp.Rational(7,20),T:sp.Rational(2,5)}
p2={H:sp.Rational(1,2),r:sp.Rational(63,250),T:sp.Rational(71,125)}
for name,(n,D) in nd.items():
    if D==0:
        expr=L*targets[name]
        check("collateral nonidentification "+name,expr.subs(p1),expr.subs(p2))
assert sp.simplify((L*targets["parent_child"]).subs(p1)-(L*targets["parent_child"]).subs(p2))!=0
l1,l2,rr,cc=sp.symbols("l1 l2 rr cc",positive=True)
check("unequal reliability attenuation",sp.sqrt(l1)*sp.sqrt(l2)*rr,sp.sqrt(l1*l2)*rr)
check("correlated error contribution",L*rr+T*cc-L*rr,T*cc)

# Closed-form covariance of the first eight people's latent phenotypes,
# derived by the existing symbolic transmission-matrix construction.
latent8=sp.zeros(8)
for i in range(8):
    for j in range(8):
        latent8[i,j]=base.red((base.q.T*base.Sigma[2*i:2*i+2,2*j:2*j+2]*base.q)[0])
latent8_fn=sp.lambdify((base.h,base.e,base.rho),latent8,"numpy")
target_fn={name:sp.lambdify((H,r),value,"numpy") for name,value in targets.items()}
def pedigree(rng,h2,rho,n):
    hv,ev=np.sqrt(h2),np.sqrt(1-h2)
    q=np.array([hv,ev])
    C=rho*np.outer(q,q)
    R=np.eye(2)+(np.sqrt(1-rho*rho)-1)*np.outer(q,q)
    sigma_s=np.sqrt((1-h2*rho)/2)
    members={"grandmother":rng.normal(size=(n,2))}
    def mate(name,parent):
        members[name]=members[parent]@C.T+rng.normal(size=(n,2))@R.T
    def child(name,p1,p2):
        z=rng.normal(size=(n,2))
        members[name]=np.column_stack(((members[p1][:,0]+members[p2][:,0])/2+sigma_s*z[:,0],z[:,1]))
    mate("grandfather","grandmother")
    child("mother","grandmother","grandfather")
    child("uncle","grandmother","grandfather")
    mate("father","mother")
    mate("aunt_by_marriage","uncle")
    child("left1","mother","father")
    child("right1","uncle","aunt_by_marriage")
    for degree in range(2,5):
        for side in ["left","right"]:
            prev=side+str(degree-1)
            mn="mate_"+prev
            mate(mn,prev)
            child(side+str(degree),prev,mn)
    names=list(members)
    X=np.column_stack([members[name] for name in names])
    P=np.column_stack([members[name]@q for name in names])
    return names,X,P
pairs={
 "spouses":("mother","father"),"parent_child":("mother","left1"),
 "full_siblings":("mother","uncle"),"grandparent":("grandmother","left1"),
 "avuncular":("uncle","left1"),"first_cousins":("left1","right1"),
 "second_cousins":("left2","right2"),"third_cousins":("left3","right3"),
 "fourth_cousins":("left4","right4"),
 "first_cousins_removed":("left2","right1"),
 "second_cousins_removed":("left3","right2"),
 "third_cousins_removed":("left4","right3")}
def simulation(h2,rho,theta,n,seed):
    rng=np.random.default_rng(seed)
    names,X,P=pedigree(rng,h2,rho,n)
    index={v:i for i,v in enumerate(names)}
    U=rng.normal(size=P.shape)
    lam=1-theta
    Y=np.sqrt(lam)*P+np.sqrt(theta)*U
    latent=np.asarray(latent8_fn(np.sqrt(h2),np.sqrt(1-h2),rho),float)
    expected=lam*latent+theta*np.eye(8)
    actual=Y[:,:8].T@Y[:,:8]/n
    se=np.sqrt((1+expected**2)/n)
    upper=np.triu_indices(8)
    max_cov_z=float(np.max(np.abs(((actual-expected)/se)[upper])))
    assert max_cov_z<6,("covariance",h2,rho,theta,max_cov_z)
    cross=U[:,:8].T@X[:,:16]/n
    max_error_cross_z=float(np.max(np.abs(cross))*np.sqrt(n))
    assert max_error_cross_z<6,("error orthogonality",max_error_cross_z)
    moments={}
    for name,(i,j) in pairs.items():
        target=lam*float(target_fn[name](h2,rho))
        estimate=float(np.mean(Y[:,index[i]]*Y[:,index[j]]))
        mcse=float(np.sqrt((1+target*target)/n))
        z=(estimate-target)/mcse
        assert abs(z)<6,(name,h2,rho,theta,z)
        moments[name]={"predicted":target,"simulated":estimate,"mc_se":mcse,"z":z}
    ychild=Y[:,index["left1"]]
    ybar=(Y[:,index["mother"]]+Y[:,index["father"]])/2
    vb=(1+lam*rho)/2
    cb=lam*h2*(1+rho)/2
    estimate_cov=float(np.mean(ychild*ybar))
    estimate_var=float(np.mean(ybar*ybar))
    estimate_child_var=float(np.mean(ychild*ychild))
    predicted_corr=cb/np.sqrt(vb)
    estimated_corr=estimate_cov/np.sqrt(estimate_var*estimate_child_var)
    corr_se=(1-predicted_corr**2)/np.sqrt(n)
    predicted_slope=cb/vb
    estimated_slope=estimate_cov/estimate_var
    slope_se=np.sqrt((1-predicted_slope**2*vb)/(n*vb))
    averages={}
    for name,pred,est,se0 in [
        ("midparent_correlation",predicted_corr,estimated_corr,corr_se),
        ("midparent_slope",predicted_slope,estimated_slope,slope_se)]:
        z=(est-pred)/se0
        assert abs(z)<6,(name,z)
        averages[name]={"predicted":float(pred),"simulated":float(est),"mc_se_asymptotic":float(se0),"z":float(z)}
    # Regression inversion checked using exact population moments, not
    # noisy logs selected for positivity.
    regression=None
    if rho>0:
        ordered=list(nd)
        design=np.array([[1,*nd[name]] for name in ordered],float)
        true_moments=np.array([lam*target_fn[name](h2,rho) for name in ordered],float)
        coef=np.linalg.lstsq(design,np.log(true_moments),rcond=None)[0]
        recovered_m=2*np.exp(coef[1])-1
        recovered_r=2*np.exp(coef[1]+coef[2])-1
        recovered_h=recovered_m/recovered_r
        recovered_theta=1-np.exp(coef[0])/recovered_h
        np.testing.assert_allclose([recovered_h,recovered_r,recovered_theta],[h2,rho,theta],atol=1e-12)
        regression={"alpha_beta_gamma":coef.tolist(),"recovered_h2_rho_noise":[recovered_h,recovered_r,recovered_theta]}
    return {"h2":h2,"rho":rho,"noise_share":theta,"n":n,"seed":seed,
      "observed_matrix_checks":36,"error_true_component_checks":128,
      "max_cov_abs_z":max_cov_z,"max_error_cross_abs_z":max_error_cross_z,
      "individual_moments":moments,"averages":averages,"population_regression_inversion":regression}
scenarios=[]
for hv in [.2,.6]:
    for rv in [.25,.75]:
        for tv in [0,.25,.6]:
            scenarios.append(simulation(hv,rv,tv,200000,19940912+len(scenarios)))
scenarios.append(simulation(.36,0,.4,200000,19940912+len(scenarios)))
# Counterexamples to a common individual-error multiplier.
rng=np.random.default_rng(19949999)
names,X,P=pedigree(rng,.36,.35,200000)
P=P[:,:8]; n=P.shape[0]
latent=np.asarray(latent8_fn(.6,.8,.35),float)
reliability=np.linspace(.35,.85,8)
U=rng.normal(size=P.shape)
Y=P*np.sqrt(reliability)+U*np.sqrt(1-reliability)
expected=np.sqrt(np.outer(reliability,reliability))*latent+np.diag(1-reliability)
actual=Y.T@Y/n
se=np.sqrt((1+expected**2)/n)
unequal_max=float(np.max(np.abs((actual-expected)/se)))
assert unequal_max<6
lam=.6; theta=.4; common=.25
U=np.sqrt(common)*rng.normal(size=(n,1))+np.sqrt(1-common)*rng.normal(size=P.shape)
Y=np.sqrt(lam)*P+np.sqrt(theta)*U
error_cov=common*np.ones((8,8))+(1-common)*np.eye(8)
expected=lam*latent+theta*error_cov
actual=Y.T@Y/n
se=np.sqrt((1+expected**2)/n)
correlated_max=float(np.max(np.abs((actual-expected)/se)))
assert correlated_max<6
wrong=lam*latent
np.fill_diagonal(wrong,1)
off=np.triu_indices(8,1)
wrong_max=float(np.max(np.abs(((actual-wrong)/se)[off])))
assert wrong_max>10
report={"symbolic_checks":len(checks),"symbolic_check_names":checks,
 "base_symbolic_checks_imported":len(base.checks),"scenarios":scenarios,
 "core_pedigrees_total":sum(s["n"] for s in scenarios),
 "individual_moment_checks":len(scenarios)*len(pairs),
 "midparent_checks":len(scenarios)*2,
 "observed_covariance_checks":len(scenarios)*36,
 "error_true_component_checks":len(scenarios)*128,
 "counterexamples":{"additional_pedigrees":n,"unequal_reliability_max_abs_z":unequal_max,
 "correlated_error_max_abs_z":correlated_max,"incorrect_common_multiplier_max_abs_z":wrong_max},
 "notes":["Six-SE diagnostics are not empirical specification tests.",
 "Midparent correlation/slope SEs are asymptotic Gaussian SEs.",
 "Regression inversion uses exact population moments; no finite-sample unbiasedness is claimed.",
 "Baseline verifier's symbol e is the review's latent loading e_P.",
 "No-error and all-noise limits are symbolic checks; the model excludes the all-noise boundary lambda=0."]}
(OUT/"measurement_error_verification.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"symbolic_checks":len(checks),"base_checks":len(base.checks),
 "core_pedigrees":report["core_pedigrees_total"],
 "individual_checks":report["individual_moment_checks"],"midparent_checks":report["midparent_checks"],
 "max_individual_abs_z":max(abs(m["z"]) for s in scenarios for m in s["individual_moments"].values()),
 "max_midparent_abs_z":max(abs(m["z"]) for s in scenarios for m in s["averages"].values()),
 "max_cov_abs_z":max(s["max_cov_abs_z"] for s in scenarios),
 "max_error_cross_abs_z":max(s["max_error_cross_abs_z"] for s in scenarios),
 "counterexamples":report["counterexamples"]},indent=2))
