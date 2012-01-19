'''
Created on 2012-01-19

@author: innovation
'''
import unittest

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
        print priors
        params = JointSnvMixParameters()        
        print params
        
        model = JointSnvMixModel(priors, params, model='jsm1')
        
        print model.params
        
        sim = JointSnvMixSimulator()
        counts, labels = sim.draw_sample(10000000)
        
        data = [JointBinaryCountData(*x) for x in counts]
        model.fit(data, verbose=True)
        
        
#        print model.params
#        q = [0] * 100000
#        r = [0] * 100000     
#        data = JointBinaryQualityData(q, r, q, r)
#        print model.predict(data)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
