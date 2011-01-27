'''
Created on 2010-07-26

@author: Andrew Roth
'''
import os
import tempfile
import unittest


from tables import openFile

from pyleup.file_formats.bpil import convert_pileup_to_bpil

class Test( unittest.TestCase ):
    def test_conversion( self ):
        pileup_file_name = "../data/ex1.pileup"
        pileup_file = open( pileup_file_name )

#        binary_pileup_file_name = "../data/ex1.bpil"
        binary_pileup_file_name = tempfile.mkstemp()[1]

        convert_pileup_to_bpil( pileup_file, binary_pileup_file_name )

        bpil_file = openFile( binary_pileup_file_name, 'r' )

        self.assertTrue( 'pileup' in bpil_file.root )
        
        self.assertTrue( 'chr1' in bpil_file.root.pileup )
        
        self.assertTrue( 'chr2' in bpil_file.root.pileup )
        
        pileup_file.close()
        bpil_file.close()
        
        os.unlink( binary_pileup_file_name )

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
