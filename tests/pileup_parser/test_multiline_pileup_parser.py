'''
Created on 2010-07-26

@author: Andrew Roth
'''
import unittest

from pyleup.pileup.parsers import PileupMultiLineParser

class Test( unittest.TestCase ):
    def setUp( self ):
        pileup_file = open( '../../../data/ex1.pileup', 'r' )

        self.lines = [x for x in pileup_file]

        self.parser = PileupMultiLineParser()

        self.chr_names = ['chr1', 'chr2']
        
        pileup_file.close()

    def test_chr_names( self ):
        columns = self.parser.parse_lines( self.lines )

        for column in columns:
            self.assertTrue( column.chr_name in self.chr_names )

    def test_numer_columns_match( self ):
        columns = self.parser.parse_lines( self.lines )
        
        self.assertEqual( len( columns ), len( self.lines ) )

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
