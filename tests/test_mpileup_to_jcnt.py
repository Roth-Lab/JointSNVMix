'''
Created on 2010-08-31

@author: Andrew Roth
'''
import unittest

from joint_snv_mix.pre_processing.mpileup_to_jcnt import get_jcnt_entry

class Test(unittest.TestCase):
    def setUp(self):
        min_depth = 1
        min_qual = 20
        
        row = {}
        row['chr_coord'] = '12345'
        row['ref_base'] = 'A'
        row['normal_depth'] = '4'
        row['normal_call_string'] = '...t'
        row['normal_base_qual_string'] = 'AAAA'
        row['tumour_depth'] = '4'
        row['tumour_call_string'] = '...c'
        row['tumour_base_qual_string'] = 'AAAA'
        
        self.jcnt_entry = get_jcnt_entry(row, min_depth, min_qual)
    
    def test_normal_matches_tumour_variant(self):        
        self.assertEqual(self.jcnt_entry[2], 'C')
        self.assertEqual(self.jcnt_entry[3], 'C')
        
    def test_normal_counts_are_correct(self):        
        self.assertEqual(self.jcnt_entry[4], 3)
        self.assertEqual(self.jcnt_entry[5], 0)
    
    def test_tumour_counts_are_correct(self):        
        self.assertEqual(self.jcnt_entry[6], 3)
        self.assertEqual(self.jcnt_entry[7], 1)
        
    def test_coord_correct(self):
        self.assertEqual(self.jcnt_entry[0], 12345)
        



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
