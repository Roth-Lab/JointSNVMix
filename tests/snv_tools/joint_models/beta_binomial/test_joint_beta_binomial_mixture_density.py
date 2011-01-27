'''
Created on 2010-09-01

@author: andrew
'''
import unittest

import numpy as np

from pyleup.snv_tools.joint_models.beta_binomial.shared_parameters import JointMixtureDensityFactory
from pyleup.snv_tools.joint_models.data import JointData

class Test( unittest.TestCase ):
    def setUp( self ):
        self.parameters = {}
        self.parameters['mix_weights'] = np.array( [1 / 9.0] * 9 )
        self.parameters['alpha'] = np.array( [9, 5, 1] )
        self.parameters['beta'] = np.array( [1, 5, 9] )
        
        factory = JointMixtureDensityFactory()        
        self.mixture_density = factory.get_mixture_density( self.parameters ) 
    
    def test_max_aa_aa_likelihood( self ):
        data = self.get_data( 50, 50, 20, 20 )
                
        component_conditional_log_likelihoods = self.mixture_density.get_component_conditional_log_likelihoods( data )
        
        self.assertEqual( np.argmax( component_conditional_log_likelihoods ), 0 )
        
    def test_max_aa_ab_likelihood( self ):
        data = self.get_data( 50, 50, 10, 20 )
                
        component_conditional_log_likelihoods = self.mixture_density.get_component_conditional_log_likelihoods( data )
        
        self.assertEqual( np.argmax( component_conditional_log_likelihoods ), 1 )
        
    def test_max_aa_bb_likelihood( self ):
        data = self.get_data( 50, 50, 1, 20 )
                
        component_conditional_log_likelihoods = self.mixture_density.get_component_conditional_log_likelihoods( data )
        
        self.assertEqual( np.argmax( component_conditional_log_likelihoods ), 2 )
        
    def test_max_bb_bb_likelihood( self ):
        data = self.get_data( 3, 50, 1, 20 )
                
        component_conditional_log_likelihoods = self.mixture_density.get_component_conditional_log_likelihoods( data )
        
        self.assertEqual( np.argmax( component_conditional_log_likelihoods ), 8 )
        
    def test_get_likelihood_returns_scalar( self ):
        data = self.get_data( 50, 50, 20, 20 )
        
        log_likelihood = self.mixture_density.get_log_likelihood( data )
        
        self.assertTrue( np.isscalar( log_likelihood ) )
                
        
    def get_data( self, k_N, n_N, k_T, n_T ):
        counts = {}
        counts['normal'] = np.array( [( k_N, n_N ),
                                      ( k_N, n_N )] )
        counts['tumour'] = np.array( [( k_T, n_T ),
                                      ( k_T, n_T )] )       
        
        data = JointData( counts )
        
        return data



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
