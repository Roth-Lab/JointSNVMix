'''
Created on 2010-08-30

@author: andrew
'''
import unittest

import numpy as np

from pyleup.snv_tools.joint_models.data import JointData


class Test( unittest.TestCase ):


    def setUp( self ):
        self.counts = {}
        self.counts['normal'] = np.array( [
                                              [0, 10],
                                              [1, 11]
                                              ] )
        
        self.counts['tumour'] = np.array( [
                                              [5, 20],
                                              [6, 21]
                                              ] )
        
        self.data = JointData( self.counts )

    def test_get_ref_matches( self ):
        self.check_ref_matches( 'normal' )
        self.check_ref_matches( 'tumour' )
        
    def test_get_non_ref_matches( self ):
        self.check_non_ref_matches( 'normal' )
        self.check_non_ref_matches( 'tumour' )
        
    def test_get_read_depth( self ):
        self.check_read_depth( "normal" )
        self.check_read_depth( "tumour" )
                
    def check_arrays_are_equal( self, array_1, array_2 ):
        return np.all( array_1 == array_2 )
    
    def check_ref_matches( self, genome ):
        ref_matches = self.data.get_ref_matches( genome )
        
        expected_ref_matches = self.counts[genome][:, 0]
        
        self.assertTrue( self.check_arrays_are_equal( ref_matches, expected_ref_matches ) )
        
    def check_non_ref_matches( self, genome ):
        non_ref_matches = self.data.get_non_ref_matches( genome )
        
        expected_non_ref_matches = self.counts[genome][:, 1] - self.counts[genome][:, 0]
        
        self.assertTrue( self.check_arrays_are_equal( non_ref_matches, expected_non_ref_matches ) )

    def check_read_depth( self, genome ):
        read_depth = self.data.get_read_depth( genome )
        
        expected_read_depth = self.counts[genome][:, 1]
        
        self.assertTrue( self.check_arrays_are_equal( read_depth, expected_read_depth ) )

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
