'''
Created on 2010-08-31

@author: Andrew Roth
'''
import unittest

import numpy as np

from joint_snv_mix.classification.utils.log_pdf import log_beta_binomial_likelihood

class Test(unittest.TestCase):
    def get_likelihood(self, k, n, alpha, beta):
        alpha = np.array(alpha)
        beta = np.array(beta)
        
        k = np.array(k)
        n = np.array(n)
        
        return log_beta_binomial_likelihood(k, n, alpha, beta)
    
    def test_likelihood_same_alpha_beta(self):
        log_likelihood = self.get_likelihood(10, 20, 5, 5)
        
        self.assertAlmostEqual(log_likelihood, -14.42888, 5)
        
    def test_likelihood_different_alpha_beta(self):
        log_likelihood = self.get_likelihood(100, 300, 1, 9)
        
        self.assertAlmostEqual(log_likelihood, -194.6557, 5)
        
    def test_likelihood_big_alpha_beta(self):
        log_likelihood = self.get_likelihood(1, 5, 9e6, 1e6)
        
        self.assertAlmostEqual(log_likelihood, -9.315696, 5)



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
