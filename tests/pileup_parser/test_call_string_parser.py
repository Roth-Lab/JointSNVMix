'''
Created on 2010-07-26

@author: Andrew Roth
'''
import unittest

from pileup.parsers import CallStringParser


class Test( unittest.TestCase ):
    def setUp( self ):
        self.parser = parser = CallStringParser()
        
    def check_ref_call_string( self, ref_base, call_string ):
        parsed_bases = self.parser.parse_call_string( ref_base, call_string )

        self.assertEqual( ref_base, parsed_bases[0] )
        
    def check_non_ref_call_string( self, ref_base, mismatch_base, call_string ):
        parsed_bases = self.parser.parse_call_string( ref_base, call_string )

        self.assertEqual( mismatch_base, parsed_bases[0] )
    
    def test_single_line_forward_strand_match( self ):        
        self.check_ref_call_string( 'A', '.' )

    def test_single_line_reverse_strand_match( self ):
        self.check_ref_call_string( 'A', ',' )

    def test_single_line_forward_strand_mismatch( self ):
        self.check_non_ref_call_string( 'A', 'C', 'C' )

    def test_single_line_reverse_strand_mismatch( self ):
        self.check_non_ref_call_string( 'A', 'C', 'c' )
  
    def test_single_line_start_of_read_ignore( self ):
        self.check_ref_call_string( 'A', '$.' )

    def test_single_line_end_of_read_ignore( self ):
        self.check_ref_call_string( 'A', '^~.' )

    def test_single_line_single_digit_insertion( self ):
        self.check_ref_call_string( 'A', "+2cc." )

    def test_single_line_multi_digit_insertion( self ):
        call_string = "+100" + "c" * 100 + "."

        self.check_ref_call_string( 'A', call_string )
        
    def test_single_line_single_digit_deletion( self ):
        self.check_ref_call_string( 'A', "-3gct." )

    def test_single_line_multi_digit_deletion( self ):        
        call_string = "-100" + "t" * 100 + "."
        
        self.check_ref_call_string( 'A', call_string )

    def test_unkown_char_exception( self ):
        self.assertRaises( Exception, self.parser.parse_call_string, ( 'A', 'Z' ) )
        
    def test_single_line_multi_digit_insertion_with_leading_char( self ):
        call_string = ".+100" + "c" * 100 + "."

        parsed_bases = self.parser.parse_call_string( 'A', call_string )
        
        self.assertEqual(parsed_bases, ['A','A'])




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
