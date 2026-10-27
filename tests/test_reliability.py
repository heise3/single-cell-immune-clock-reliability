"""Synthetic independent checks; no participant files or empirical output writes."""
from pathlib import Path
import importlib.util, unittest
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('reliability',P/'outputs/round2_A1/compute_reliability.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

class ReliabilityChecks(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(20261025)
        self.g=pd.DataFrame([(f'd{i}',float(25+i*4),j) for i in range(12) for j in range(2+i%3)],columns=['donor','age','visit'])
        self.y=4+.7*self.g.age.to_numpy()+rng.normal(0,2,len(self.g))
    def oracle(self,mult):
        g=self.g;ids=sorted(g.donor.unique());out=np.empty(len(g))
        for held in ids:
            indices=[];weights=[]
            for d,ncopy in zip(ids,mult):
                ix=np.flatnonzero(g.donor.eq(d))
                if d!=held:
                    indices.extend(np.tile(ix,int(ncopy)));weights.extend([1/len(ix)]*(len(ix)*int(ncopy)))
            indices=np.asarray(indices,int);sw=np.sqrt(weights)
            X=np.c_[np.ones(len(indices)),g.age.to_numpy()[indices]]
            coef=np.linalg.lstsq(X*sw[:,None],self.y[indices]*sw,rcond=None)[0]
            at=np.flatnonzero(g.donor.eq(held));out[at]=self.y[at]-np.c_[np.ones(len(at)),g.age.to_numpy()[at]]@coef
        return out
    def test_absolute_agreement_hand_example(self):
        a=np.array([[1,2],[2,3],[3,4],[4,5]],float)
        icc,sem,mdc=m.icc_sem(a)
        self.assertAlmostEqual(icc,10/13)
        self.assertAlmostEqual(sem,np.std(a.ravel(),ddof=1)*np.sqrt(3/13))
        self.assertAlmostEqual(mdc,1.96*np.sqrt(2)*sem)
    def test_perfect_and_undefined(self):
        np.testing.assert_allclose(m.icc_sem(np.c_[np.arange(8),np.arange(8)]),[1,0,0],atol=1e-12)
        self.assertTrue(np.isnan(m.icc_sem(np.ones((8,2)))).all())
    def test_constant_common_shift(self):
        a=self.y[:24].reshape(-1,2)
        np.testing.assert_allclose(m.icc_sem(a),m.icc_sem(a+70),rtol=1e-12,atol=1e-12)
    def test_lodo_unequal_visits(self):
        np.testing.assert_allclose(m.lodo(self.g,self.y),self.oracle(np.ones(12,int)),atol=1e-10)
    def test_all_bootstrap_copies_excluded(self):
        b=np.random.default_rng(20261025).multinomial(12,np.ones(12)/12,size=20)
        actual=m.batched_lodo(self.g,self.y,b)
        for i,mult in enumerate(b):np.testing.assert_allclose(actual[i],self.oracle(mult),atol=1e-9)
    def test_held_out_response_does_not_fit_itself(self):
        ix=self.g.donor.eq('d3').to_numpy();r1=m.lodo(self.g,self.y);y2=self.y.copy();y2[ix]+=100
        r2=m.lodo(self.g,y2);np.testing.assert_allclose(r2[ix]-r1[ix],100,atol=1e-10)

if __name__=='__main__':unittest.main()
