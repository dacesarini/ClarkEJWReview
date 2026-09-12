"""Verify the Fisher-model review section. Run from any working directory."""
from pathlib import Path
import json
import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parents[1] / "_build"
OUT.mkdir(exist_ok=True)
h, e, rho = sp.symbols("h e rho", real=True)
q = sp.Matrix([h, e])
I = sp.eye(2)
checks = []
def red(x):
    return sp.rem(sp.Poly(sp.expand(x), e), sp.Poly(e**2+h**2-1, e)).as_expr()
def matred(m):
    return m.applyfunc(red)
def check(label, actual, expected):
    if isinstance(actual, sp.MatrixBase):
        residual = matred(actual-expected)
        assert residual == sp.zeros(*residual.shape), (label, residual)
    else:
        residual = red(actual-expected)
        assert residual == 0, (label, residual)
    checks.append(label)

C = rho*q*q.T
k = (1+rho*h*h)/2
d = rho*h*e/2
seg = (1-rho*h*h)/2
check("phenotype variance", (q.T*q)[0], 1)
check("spouse phenotype covariance", (q.T*C*q)[0], rho)
check("spouse genetic covariance", C[0,0], rho*h*h)
check("PH Gaussian conditional cross covariance given wife's phenotype",
      C-q*(q.T*C), sp.zeros(2))
check("PH Gaussian conditional cross covariance given husband's phenotype",
      C-(C*q)*q.T, sp.zeros(2))
check("mate residual covariance", I-C*C, I-rho*rho*q*q.T)
check("determinant of full spouse covariance",
      sp.det(I.row_join(C).col_join(C.T.row_join(I))), 1-rho*rho)
Sigma = I.row_join(C).col_join(C.T.row_join(I))
people = {"grandmother":0, "grandfather":2}
def append_pair(L, D):
    global Sigma
    cross = matred(Sigma*L.T)
    block = matred(L*Sigma*L.T+D)
    start = Sigma.rows
    Sigma = Sigma.row_join(cross).col_join(cross.T.row_join(block))
    return start
def child(name, p1, p2):
    L = sp.zeros(2, Sigma.rows)
    L[0,people[p1]] = sp.Rational(1,2)
    L[0,people[p2]] = sp.Rational(1,2)
    people[name] = append_pair(L, sp.diag(seg, 1))
def mate(name, parent):
    L = sp.zeros(2, Sigma.rows)
    j = people[parent]
    L[:,j:j+2] = C
    people[name] = append_pair(L, matred(I-C*C))
child("mother","grandmother","grandfather")
child("uncle","grandmother","grandfather")
mate("father","mother")
mate("aunt_by_marriage","uncle")
child("child","mother","father")
child("cousin","uncle","aunt_by_marriage")
def block(a,b):
    i,j=people[a],people[b]
    return Sigma[i:i+2,j:j+2]
def pcov(a,b):
    return (q.T*block(a,b)*q)[0]
for name in people:
    check(name+" marginal component covariance",block(name,name),I)
    check(name+" marginal phenotype variance",pcov(name,name),1)
for a,b in [("grandmother","grandfather"),("mother","father"),("uncle","aunt_by_marriage")]:
    check(a+"/"+b+" spouse block",block(a,b),C)
check("stationarity", k+seg,1)
check("full sibling component block",block("mother","uncle"),sp.diag(k,0))
check("parent child component block",block("mother","child"),sp.Matrix([[k,0],[d,0]]))
check("grandparent child component block",block("grandmother","child"),sp.Matrix([[k*k,0],[d*k,0]]))
check("avuncular component block",block("uncle","child"),sp.diag(k*k,0))
check("first cousin component block",block("child","cousin"),sp.diag(k**3,0))
targets = {
    "spouses":("mother","father",rho),
    "parent_child":("mother","child",h*h*(1+rho)/2),
    "full_siblings":("mother","uncle",h*h*k),
    "grandparent":("grandmother","child",h*h*(1+rho)*k/2),
    "avuncular":("uncle","child",h*h*k*k),
    "first_cousins":("child","cousin",h*h*k**3),
}
for name,(a,b,target) in targets.items():
    check(name+" phenotype correlation",pcov(a,b),target)
check("parent minus sibling gap",
      pcov("mother","child")-pcov("mother","uncle"),
      h*h*rho*(1-h*h)/2)
for name,(_,_,target) in targets.items():
    random_expected={"spouses":0,"parent_child":h*h/2,
                     "full_siblings":h*h/2,"grandparent":h*h/4,
                     "avuncular":h*h/4,"first_cousins":h*h/8}[name]
    check(name+" random mating limit",target.subs(rho,0),random_expected)

# Law-of-total-covariance check: extra mate residual association would
# add h^2 * delta/4 to the cousin phenotype covariance.
delta=sp.symbols("delta",real=True)
check("cousin sensitivity to additional mate genetic covariance",
      h*h*(k**3+delta/4)-h*h*k**3,h*h*delta/4)

numeric_cov = sp.lambdify((h,e,rho),Sigma,"numpy")
numeric_targets = {name:sp.lambdify((h,e,rho),v[2],"numpy")
                   for name,v in targets.items()}
def simulation(h2, r, n, seed):
    rng=np.random.default_rng(seed)
    hv,ev=np.sqrt(h2),np.sqrt(1-h2)
    qv=np.array([hv,ev])
    Cv=r*np.outer(qv,qv)
    # Symmetric square root of I-C*C: eigenvalues 1 and sqrt(1-r^2).
    R=np.eye(2)+(np.sqrt(1-r*r)-1)*np.outer(qv,qv)
    segregation_sd=np.sqrt((1-r*h2)/2)
    members={}
    members["grandmother"]=rng.normal(size=(n,2))
    def draw_mate(name,parent):
        members[name]=members[parent]@Cv.T+rng.normal(size=(n,2))@R.T
    def draw_child(name,p1,p2):
        noise=rng.normal(size=(n,2))
        members[name]=np.column_stack((
            (members[p1][:,0]+members[p2][:,0])/2+segregation_sd*noise[:,0],
            noise[:,1]))
    draw_mate("grandfather","grandmother")
    draw_child("mother","grandmother","grandfather")
    draw_child("uncle","grandmother","grandfather")
    draw_mate("father","mother")
    draw_mate("aunt_by_marriage","uncle")
    draw_child("child","mother","father")
    draw_child("cousin","uncle","aunt_by_marriage")
    X=np.column_stack([members[name] for name in people])
    expected=np.asarray(numeric_cov(hv,ev,r),dtype=float)
    min_eigenvalue=float(np.linalg.eigvalsh(expected).min())
    assert min_eigenvalue > 0
    actual=X.T@X/n
    # For centered jointly normal U,V: Var(UV)=Var(U)Var(V)+Cov(U,V)^2.
    se=np.sqrt((np.outer(np.diag(expected),np.diag(expected))+expected**2)/n)
    standardized=(actual-expected)/se
    upper=np.triu_indices(expected.shape[0])
    component_max=float(np.max(np.abs(standardized[upper])))
    assert component_max < 6, ("component Monte Carlo check",h2,r,component_max)
    moments={}
    for name,(a,b,_) in targets.items():
        pa=members[a]@qv
        pb=members[b]@qv
        estimate=float(np.mean(pa*pb))
        target=float(numeric_targets[name](hv,ev,r))
        mcse=float(np.sqrt((1+target*target)/n))
        z=(estimate-target)/mcse
        assert abs(z)<6,(name,h2,r,z)
        moments[name]={"predicted":target,"simulated":estimate,"mc_se":mcse,"z":z}
    return {"h2":h2,"rho":r,"families":n,"seed":seed,
            "minimum_eigenvalue":min_eigenvalue,
            "component_covariance_checks":len(upper[0]),
            "component_max_abs_z":component_max,"phenotype_moments":moments}
scenarios=[]
for h2 in [.1,.36,.8]:
    for r in [0,.35,.8]:
        scenarios.append(simulation(h2,r,250000,20260912+len(scenarios)))
report={
 "symbolic_checks":len(checks),"symbolic_check_names":checks,
 "pedigree_members_in_matrix_order":list(people),
 "simulation_scenarios":scenarios,
 "families_total":sum(s["families"] for s in scenarios),
 "component_covariance_checks_total":sum(s["component_covariance_checks"] for s in scenarios),
 "phenotype_checks_total":len(scenarios)*len(targets),
 "diagnostic":"Six-standard-error threshold, with known zero means and Gaussian product variances; not a hypothesis test or proof.",
 "scope":"Explicit stationary Gaussian moment construction; not a simulation of allele frequencies or convergence of a multilocus population."
}
(OUT/"fisher_model_verification.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
chosen=next(s for s in scenarios if s["h2"]==.36 and s["rho"]==.35)
print(json.dumps({"symbolic_checks":len(checks),"families_total":report["families_total"],
 "component_checks":report["component_covariance_checks_total"],
 "phenotype_checks":report["phenotype_checks_total"],
 "largest_component_abs_z":max(s["component_max_abs_z"] for s in scenarios),
 "largest_phenotype_abs_z":max(abs(m["z"]) for s in scenarios for m in s["phenotype_moments"].values()),
 "example":chosen},indent=2))
