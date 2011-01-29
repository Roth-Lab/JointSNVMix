'''
Created on 2010-07-26

@author: Andrew Roth
'''
import unittest

from pileup.parsers import PileupLineParser

class Test( unittest.TestCase ):
    def setUp( self ):
        self.chr_name = "chr1"

        self.position = 1234

        self.ref_base = "T"

        self.read_depth = 2

        self.call_string = "^~.^~g"
        
        self.bases = ['T', 'G']

        self.base_qual_string = "B@"
        
        self.base_quals = [33, 31]

        self.map_qual_string = "~:"
        
        self.map_quals = [93, 25]

        self.pileup_line_parts = [
                                  self.chr_name,
                                  str( self.position ),
                                  self.ref_base,
                                  str( self.read_depth ),
                                  self.call_string,
                                  self.base_qual_string,
                                  self.map_qual_string
                                  ]

        self.pileup_line = ""

        for part in self.pileup_line_parts:
            self.pileup_line += part + "\t"

        # Switch trailing tab to newline.
        self.pileup_line = self.pileup_line[:-1] + "\n"

    def test_chr_name_parsing( self ):
        parser = PileupLineParser()

        parsed_chr_name = parser.parse_line( self.pileup_line ).chr_name

        self.assertEqual( self.chr_name, parsed_chr_name )
        
    def test_position_parsing( self ):
        parser = PileupLineParser()

        parsed_position = parser.parse_line( self.pileup_line ).position

        self.assertEqual( self.position, parsed_position )

    def test_ref_base_parsing( self ):
        parser = PileupLineParser()

        parsed_ref_base = parser.parse_line( self.pileup_line ).ref_base

        self.assertEqual( self.ref_base, parsed_ref_base )
        
    def test_bases_parsing( self ):
        parser = PileupLineParser()

        parsed_bases = parser.parse_line( self.pileup_line ).bases

        self.assertEqual( self.bases, parsed_bases )
    
    def test_base_quals_parsing( self ):
        parser = PileupLineParser()

        parsed_base_quals = parser.parse_line( self.pileup_line ).base_quals

        self.assertEqual( self.base_quals, parsed_base_quals )
        
    def test_map_quals_parsing( self ):
        parser = PileupLineParser()

        parsed_map_quals = parser.parse_line( self.pileup_line ).map_quals

        self.assertEqual( self.map_quals, parsed_map_quals )
        
    def test_problem_string( self ):
        string = "chr22   16050581        a       1       ^F.     G       F"
        string = string.split()
        string = "\t".join(string)
        
        parser = PileupLineParser()
        
        result = parser.parse_line( string )
        
        print result.map_quals

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
