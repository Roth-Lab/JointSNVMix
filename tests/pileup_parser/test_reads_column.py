'''
Created on 2010-07-26

@author: Andrew Roth
'''

import unittest
from tests.pileup.generators import ReadsColumnGenerator

class Test( unittest.TestCase ):
    def setUp( self ):
        self.generator = ReadsColumnGenerator()

    def test_get_counts( self ):
        expected_values, column = self.generator.create_read( 1000 )

        expect = self.generator.get_counts( expected_values['bases'] )

        result = column.get_counts()

        self.assertEqual( expect, result )

    def test_get_counts_base_qual_threshold( self ):
        expected_values, column = self.generator.create_read( 1000 )

        expect = self.generator.get_counts( expected_values['bases'], min_base_qual = 10 )

        result = column.get_counts( min_base_qual = 10 )

        self.assertEqual( expect, result )
        
    def test_get_counts_map_qual_threshold( self ):
        expected_values, column = self.generator.create_read( 1000 )

        expect = self.generator.get_counts( expected_values['bases'], min_map_qual = 10 )

        result = column.get_counts( min_map_qual = 10 )

        self.assertEqual( expect, result )
        
    def test_get_counts_both_qual_threshold( self ):
        expected_values, column = self.generator.create_read( 1000 )

        expect = self.generator.get_counts( expected_values['bases'],  min_base_qual = 20, min_map_qual = 10 )

        result = column.get_counts(  min_base_qual = 20, min_map_qual = 10 )

        self.assertEqual( expect, result )




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
