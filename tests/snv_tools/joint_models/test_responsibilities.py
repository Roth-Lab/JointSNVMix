'''
Created on 2010-08-31

@author: bcca
'''
import unittest

import numpy as np

import pyleup.constants as constants

from pyleup.snv_tools.joint_models.responsibilities import Responsibilities

class Test( unittest.TestCase ):
    def setUp( self ):
        self.nrows = 100
        self.ncols = 9

        self.shape = ( self.nrows, self.ncols )

        self.component_names = constants.joint_genotypes

        self.responsibilities = Responsibilities( self.nrows, self.component_names )

        self.component_log_likelihoods = np.random.random( self.shape )

        self.responsibilities.update( self.component_log_likelihoods )

    def test_normalises_across_columns( self ):
        resp_array = self.responsibilities.to_array()

        sum_across_columns = resp_array.sum( axis = 1 )

        is_close_to_one = np.allclose( sum_across_columns, 1, 0, 1e-12 )

        self.assertTrue( is_close_to_one )

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
