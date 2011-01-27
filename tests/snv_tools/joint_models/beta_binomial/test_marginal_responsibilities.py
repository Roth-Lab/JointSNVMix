'''
Created on 2010-09-01

@author: andrew
'''
import unittest

import numpy as np

import pyleup.constants as constants
from pyleup.snv_tools.joint_models.responsibilities import Responsibilities
from pyleup.snv_tools.joint_models.beta_binomial.shared_parameters import MarginalResponsibilities


class Test( unittest.TestCase ):
    def setUp( self ):
        self.nrows = 100
        self.ncols = 9

        self.shape = ( self.nrows, self.ncols )

        self.component_names = constants.joint_genotypes

        self.responsibilities = Responsibilities( self.nrows, self.component_names )

        self.component_log_likelihoods = np.random.random( self.shape )

        self.responsibilities.update( self.component_log_likelihoods )
        
        self.marginal_resp = MarginalResponsibilities( self.responsibilities )  
            
    def test_marginal_responsibilities_shape( self ):
        self.assertTrue( self.marginal_resp.to_array( 'normal' ).shape == ( self.nrows, 3 ) )
        self.assertTrue( self.marginal_resp.to_array( 'tumour' ).shape == ( self.nrows, 3 ) )

    def test_normal_marginalisation( self ):
        normal_indices = [( 0, 1, 2 ), ( 3, 4, 5 ), ( 6, 7, 8 )]

        expected_marginals = self.get_expected_marginals( normal_indices )

        self.assertTrue( np.allclose( self.marginal_resp.to_array( 'normal' ), expected_marginals, 0, 1e-10 ) )

    def test_tumour_marginalisation( self ):
        tumour_indices = [ ( 0, 3, 6 ), ( 1, 4, 7 ), ( 2, 5, 8 )]
        
        expected_marginals = self.get_expected_marginals( tumour_indices )

        self.assertTrue( np.allclose( self.marginal_resp.to_array( 'tumour' ), expected_marginals, 0, 1e-10 ) )

    def get_expected_marginals( self, indices ):
        resp = self.responsibilities.to_array()

        marginal_columns = []

        for index_row in indices:
            marginal_columns.append( resp[:, index_row].sum( axis = 1 ) )

        expected_marginals = np.column_stack( marginal_columns )

        return expected_marginals


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
