'''
Created on 2010-09-01

@author: andrew
'''
import unittest

import numpy as np

from pyleup.snv_tools.joint_models.beta_binomial.densities import BetaBinomialDensity, \
    JointBetaBinomialDensity
    
from pyleup.snv_tools.joint_models.data import JointData

class Test( unittest.TestCase ):
    def test_joint_likelihood_same_density( self ):
        '''
        In R
        >>> library(VGAM)
        >>> x = dbetabin.ab(30,50,1,9) / choose(50,30)
        >>> y = dbetabin.ab(10,20,1,9) / choose(20,10)
        >>> log_likelihood = log(x) + log(y)
        '''
         
        normal_density = BetaBinomialDensity( 1, 9 )        
        tumour_density = BetaBinomialDensity( 1, 9 )
        
        joint_density = JointBetaBinomialDensity( normal_density, tumour_density )
        
        k_N = 30
        n_N = 50
        
        k_T = 10
        n_T = 20
        
        counts = {}
        counts['normal'] = np.atleast_2d( np.array( [k_N, n_N] ) )
        counts['tumour'] = np.atleast_2d( np.array( [k_T, n_T] ) )       
        
        data = JointData( counts )
        
        log_likelihood = joint_density.get_log_likelihood( data )
        
        self.assertAlmostEqual( log_likelihood, -57.34858, 5 )
        
    def test_joint_likelihood_different_densities( self ):
        normal_density = BetaBinomialDensity( 5, 5 )        
        tumour_density = BetaBinomialDensity( 1, 9 )
        
        joint_density = JointBetaBinomialDensity( normal_density, tumour_density )
        
        k_N = 30
        n_N = 50
        
        k_T = 10
        n_T = 20
        
        counts = {}
        counts['normal'] = np.atleast_2d( np.array( [k_N, n_N] ) )
        counts['tumour'] = np.atleast_2d( np.array( [k_T, n_T] ) )       
        
        data = JointData( counts )
        
        log_likelihood = joint_density.get_log_likelihood( data )
        
        self.assertAlmostEqual( log_likelihood, -52.28251, 5 )
    
    def test_joint_pdf_different_densities( self ):
        normal_density = BetaBinomialDensity( 5, 5 )        
        tumour_density = BetaBinomialDensity( 1, 9 )
        
        joint_density = JointBetaBinomialDensity( normal_density, tumour_density )
        
        k_N = 30
        n_N = 50
        
        k_T = 10
        n_T = 20
        
        counts = {}
        counts['normal'] = np.atleast_2d( np.array( [k_N, n_N] ) )
        counts['tumour'] = np.atleast_2d( np.array( [k_T, n_T] ) )       
        
        data = JointData( counts )
        
        log_pdf = joint_density.get_log_pdf( data )
        
        self.assertAlmostEqual( log_pdf, -8.671803, 5 )
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
