'''
Created on 2012-01-19

@author: innovation
'''
import unittest

import numpy as np
from tests.simualtors.joint_binomial import JointSnvMixSimulator

from joint_snv_mix.counter import JointBinaryCountData, JointBinaryQualityData
from joint_snv_mix.models.joint_snv_mix import JointSnvMixModel, JointSnvMixPriors, JointSnvMixParameters

class Test(unittest.TestCase):
    def test_init(self):
        mu = (
              {'alpha' : 100, 'beta' : 2},
              {'alpha' : 50, 'beta' : 50},
              {'alpha' : 2, 'beta' : 100}
              )
        
        priors = JointSnvMixPriors(mu_N=mu, mu_T=mu)
        params = JointSnvMixParameters()        
        
        model = JointSnvMixModel(priors, params, model='jsm2')
        
#        sim = JointSnvMixSimulator(mu_T=[0.9, 0.3, 0.01])
#        counts, labels = sim.draw_sample(100000)
#        
#        data = [JointBinaryCountData(*x) for x in counts]
#        model.fit(data, verbose=True)
#        
#        correct = 0
#        
#        for x, t in zip(data, labels):
#            p = model.predict(x)
#            
#            l = p.index(max(p))
#            t = np.argmax(t)
#            
#            if t == l:
#                correct += 1
#        
#        print correct
#        
        
        print model.params
        q = [0] * 100
        r = [1] * 100
        data = [JointBinaryQualityData(q, r, q, r) for _ in range(100000)]
        model.fit(data, verbose=True)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
