'''
Created on 2011-03-10

@author: andrew
'''
import unittest
from joint_snv_mix.pre_processing.bam_to_jcnt import convert_positions_to_regions


class Test(unittest.TestCase):
    def setUp(self):
        test_file_name = 'data/test.pos'
        
        self.regions = convert_positions_to_regions(test_file_name)
        
    def test_first_region(self):
        self.assertEqual(['1', 1872, 1919], self.regions[0])

    def test_last_region(self):        
        self.assertEqual(['Y', 57410923, 57411022], self.regions[-1])
        
    def test_other_regions(self):        
        self.assertEqual(['1', 2041, 2089], self.regions[1])
        self.assertEqual(['1', 2475, 2477], self.regions[2])
        self.assertEqual(['X', 110409, 110508], self.regions[3])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
