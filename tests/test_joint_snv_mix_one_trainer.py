'''
Created on 2011-09-07

@author: andrew
'''
import unittest

import random

from joint_snv_mix.trainers.joint_snv_mix import JointSnvMixOneData, JointSnvMixParameters, JointSnvMixPriors, \
                                                 JointSnvMixOneModel, JointSnvMixOneCpt
                                                 
from tests.simualtors.joint_binomial import JointSnvMixSimulator                                                 

class Test(unittest.TestCase):
    def test_data(self):
        normal_counts = (100, 100)
        tumour_counts = (1000, 500)
        
        data = JointSnvMixOneData(normal_counts, tumour_counts)
        
        self.assertSequenceEqual(normal_counts, data.normal)
        self.assertSequenceEqual(tumour_counts, data.tumour)
    
    def test_test_init_priors(self):
        priors = JointSnvMixPriors()
        
        default_mu = ((100, 2), (50, 50), (2, 100))
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)

        self.assertSequenceEqual(priors.mu_N, default_mu)
        self.assertSequenceEqual(priors.mu_T, default_mu)
        self.assertSequenceEqual(priors.pi, default_pi)
        
        print priors
    
    def test_test_init_parameters(self):
        params = JointSnvMixParameters()
        
        default_mu = (0.99, 0.5, 0.01)
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        default_pi = [x / sum(default_pi) for x in default_pi]
        
        self.assertSequenceEqual(params.mu_N, default_mu)
        self.assertSequenceEqual(params.mu_T, default_mu)
        self.assertSequenceEqual(params.pi, default_pi)
        
        print params
        
    def test_model_fit(self):
        dataset = self._get_test_dataset(1000000)
        
        params = JointSnvMixParameters()
        model = JointSnvMixOneModel(params)
        
        model.fit(dataset)
#        
    def test_cpt(self):
        normal_counts = (10, 0)
        tumour_counts = (10, 10)
        
        parms = JointSnvMixParameters()
        data = JointSnvMixOneData(normal_counts, tumour_counts)
        
        cpt = JointSnvMixOneCpt(data, parms)
        
        print cpt.resp
        print cpt.marginal
    
    def _get_test_dataset(self, size):       
        dataset = []
        
        simulator = JointSnvMixSimulator(mu_T=(0.99, 0.7, 0.001), normal_depth=100, tumour_depth=100)
        sample, labels = simulator.draw_sample(size)        
        
        for row in sample:
            normal_counts = (row[0], row[1])
            tumour_counts = (row[2], row[3])
            
            data = JointSnvMixOneData(normal_counts, tumour_counts)
            
            dataset.append(data)
        
        return dataset

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
