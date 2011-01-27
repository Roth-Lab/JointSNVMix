'''
Created on 2010-07-29

@author: bcca
'''
import unittest

from pyleup.pileup.binary_pileup import BinaryPileupReader

class Test( unittest.TestCase ):
    def setUp( self ):
        self.reader = BinaryPileupReader( '../../../data/ex1.bpil' )

    def tearDown( self ):
        self.reader.close()

    def test_get_chromosomes( self ):
        expected_chromosomes = set( ['chr1', 'chr2'] )

        self.assertEqual( expected_chromosomes, self.reader.get_chromosomes() )

    def test_get_positions( self ):
        self.assertEqual( len( self.reader.get_positions( 'chr1' ) ), 1470 )

    def test_get_column( self ):
        column = self.reader.get_column( 'chr1', 1555 )

        self.assertEqual( column.ref_base, 'C' )
        self.assertEqual( len( column.reads ), 18 )

    def test_get_column_iter( self ):
        chr_name = 'chr1'

        positions = self.reader.get_positions( chr_name )

        for ( i, column ) in enumerate( self.reader.get_column_iterator( chr_name ) ):
            self.assertEqual( column.position, positions[i] )

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
