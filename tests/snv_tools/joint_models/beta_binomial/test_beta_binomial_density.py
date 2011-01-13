'''
Created on 2010-08-31

@author: Andrew Roth
'''
import unittest

from scipy.special import gammaln

from pyleup.utils.log_pdf import log_binomial_coefficient

from pyleup.snv_tools.joint_models.beta_binomial.densities import BetaBinomialDensity

class Test( unittest.TestCase ):
    def test_likelihood_same_alpha_beta( self ):
        alpha = 5
        beta = 5
        
        k = 10
        n = 20
        
        density = BetaBinomialDensity( alpha, beta )
        
        log_likelihood = density.get_log_likelihood( k, n )
        
        self.assertAlmostEqual( log_likelihood, -14.42888, 5 )
        
    def test_likelihood_different_alpha_beta( self ):
        alpha = 1
        beta = 9
        
        k = 100
        n = 300
        
        density = BetaBinomialDensity( alpha, beta )
        
        log_likelihood = density.get_log_likelihood( k, n )
        
        self.assertAlmostEqual( log_likelihood, -194.6557, 5 )
        
    def test_likelihood_big_alpha_beta( self ):
        alpha = 9e6
        beta = 1e6
        
        k = 1
        n = 5
        
        density = BetaBinomialDensity( alpha, beta )
        
        log_likelihood = density.get_log_likelihood( k, n )
        
        self.assertAlmostEqual( log_likelihood, -9.315696, 5 )
        
    def test_pdf( self ):
        alpha = 5
        beta = 5
        
        k = 10
        n = 20
        
        density = BetaBinomialDensity( alpha, beta )
        
        log_pdf = density.get_log_pdf( k, n )
        
        self.assertAlmostEqual( log_pdf, -2.302085, 5 )



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
