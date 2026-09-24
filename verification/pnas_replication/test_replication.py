"""Numerical/unit tests; these do not certify exact replication of Clark."""
import unittest
import numpy as np
from scipy.optimize import differential_evolution
from build_inputs import binary_stats,DISTANCE,LINEAL
from estimate import fit,published,weights,log_prediction

class ReplicationTests(unittest.TestCase):
    def test_chat_ols(self):
        targets={'1780-1859':(.577604823393,.526418977812,.794448213706),
                 '1860-1919':(.576130957842,.428345146747,.253552923745)}
        for c,target in targets.items():
            r=fit(published(c));np.testing.assert_allclose([r[k] for k in ('h2','m','r')],target,atol=1e-10)

    def test_chat_wls(self):
        r=fit(published('1780-1859'),'normal_log')
        np.testing.assert_allclose([r[k] for k in ('h2','m','r')],[.621531518,.555315522,.705314004],atol=1e-8)

    def test_ph_constraint_and_nested_objective(self):
        for c in ('1780-1859','1860-1919'):
            for w in ('equal','normal_corr','normal_log'):
                u=fit(published(c),w);r=fit(published(c),w,'ph_fixed')
                self.assertAlmostEqual(r['m'],r['h2']*r['r'],12)
                self.assertGreaterEqual(r['objective_normalized']+1e-10,u['objective_normalized'])

    def test_global_optimizer_check(self):
        data=published('1860-1919');rho=np.array(data['rho']);n=np.array(data['n']);d=np.array(data['d_lin'])
        for kind in ('equal','normal_log'):
            w=weights(kind,rho,np.array(data['N']));w/=w.mean()
            def objective(p):return np.sum(w*(log_prediction(p,n,d,'ph_fixed')-np.log(rho))**2)
            global_fit=differential_evolution(objective,[(1e-6,1),(1e-6,1)],seed=2026,tol=1e-10)
            local=fit(data,kind,'ph_fixed')
            self.assertAlmostEqual(global_fit.fun,local['objective_normalized'],8)

    def test_synthetic_recovery(self):
        p=np.array([.65,.75]);n=np.array(DISTANCE);d=np.array(LINEAL)
        data=dict(rho=np.exp(log_prediction(p,n,d,'ph_fixed')),N=np.full(11,1000),n=n,d_lin=d)
        r=fit(data,'normal_log','ph_fixed');np.testing.assert_allclose([r['h2'],r['r']],p,atol=1e-8)
        it=fit(data,'normal_log','ph_fixed',iterate=True);self.assertTrue(it['converged'])

    def test_binary_se_influence_against_multinomial_delta(self):
        cells=np.array([800,70,80,50]);p=cells/cells.sum();r=binary_stats(cells)
        def corr(q):
            px=q[2]+q[3];py=q[1]+q[3]
            return (q[3]-px*py)/np.sqrt(px*(1-px)*py*(1-py))
        eps=1e-6;grad=np.array([(corr(p+eps*np.eye(4)[i])-corr(p-eps*np.eye(4)[i]))/(2*eps) for i in range(4)])
        cov=(np.diag(p)-np.outer(p,p))/cells.sum()
        self.assertAlmostEqual(r['se_binary']**2,float(grad@cov@grad),10)

    def test_scale_invariance_of_supplied_weights(self):
        data=published('1780-1859');data['weight']=np.arange(1,12)
        a=fit(data,'supplied');data['weight']=data['weight']*100;b=fit(data,'supplied')
        np.testing.assert_allclose([a[k] for k in ('a','beta','c')],[b[k] for k in ('a','beta','c')],atol=1e-12)

    def test_reject_invalid_inputs(self):
        data=published('1780-1859');data['rho']=list(data['rho']);data['rho'][0]=0
        with self.assertRaises(ValueError):fit(data)
        with self.assertRaises(ValueError):fit(published('1780-1859'),'binary_corr')

if __name__=='__main__':unittest.main(verbosity=2)
