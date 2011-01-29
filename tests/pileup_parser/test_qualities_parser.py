'''
Created on 2010-07-26

@author: Andrew Roth
'''
import unittest

from pyleup.pileup.parsers import QualityStringParser

class Test( unittest.TestCase ):
    def test_quality_values( self ):
        '''
        Check the quality parser gives the right output.
        
        Use ascii characters from 33 - (93 + 33). Corresponds to qualities from
        0 - 93
        '''
        quality_string = ""

        for i in range( 93 ):
            quality_string += chr( i + 33 )

        parser = QualityStringParser()
        
        qualities = parser.parse_quality_string( quality_string )

        for ( i, quality ) in enumerate( qualities ):
            self.assertEqual( i, quality )




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
