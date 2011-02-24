'''
Created on 2010-07-26

@author: Andrew Roth
'''
import unittest

from joint_snv_mix.pre_processing.mpileup_to_jcnt import get_bases

class Test(unittest.TestCase):
    def test_quality_values(self):
        '''
        Check the quality parser gives the right output.
        
        Use ascii characters from 33 - (93 + 33). Corresponds to qualities from
        0 - 93
        '''
        call_string = "..t."
        quality_string_int = [30, 20, 20, 1]
        qual_string = ""

        for qual_int in quality_string_int:
            qual_string += chr(qual_int + 33)
            
        bases = get_bases('A', call_string, qual_string, 30)        
        self.assertEqual(bases.tolist(), ['A'])
        
        bases = get_bases('A', call_string, qual_string, 20)        
        self.assertEqual(bases.tolist(), ['A', 'A', 'T'])
        
        bases = get_bases('A', call_string, qual_string, 1)        
        self.assertEqual(bases.tolist(), ['A', 'A', 'T', 'A'])





if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
