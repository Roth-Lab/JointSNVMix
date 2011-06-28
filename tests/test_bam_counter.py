'''
Created on 2011-06-21

@author: andrew
'''
import unittest

import os.path

import pysam

from joint_snv_mix.pre_processing.base_counter import BaseCounter

class Test(unittest.TestCase):
    def setUp(self):
        self.bam_folder = "data/bam"
    
    def get_counter(self, bam_file_base_name, min_base_qual=10, min_map_qual=10):
        bam_file_name = os.path.join(self.bam_folder, bam_file_base_name)
        
        bam_file = pysam.Samfile(bam_file_name, 'rb')        
        reader = BaseCounter(bam_file, min_base_qual, min_map_qual)

        return reader
    
    def test_refs(self):
        reader = self.get_counter('plain.bam')
        
        self.assertTrue('1' in reader.refs)
        
    def test_iter_row(self):
        reader = self.get_counter('plain.bam')
        
        iter = reader.iter_ref('1')
        
        row_1 = iter.next()        
        self.assertEqual(row_1.position, 33)
        self.assertEqual(row_1.counts['A'], 1)
        self.assertEqual(sum(row_1.counts.values()), 1)
    
    def test_deletion_correct(self):
        reader = self.get_counter('deletion.bam')
        
        for row in reader.iter_ref('1'):
            if row.position == 34:
                self.assertEqual(sum(row.counts.values()), 0)
                
            if row.position == 35:
                self.assertEqual(row.counts['C'], 1)
                self.assertEqual(row.counts['G'], 2)

    def test_insertion_correct(self):
        reader = self.get_counter('insertion.bam')
        
        for row in reader.iter_ref('1'):
            if row.position == 34:
                self.assertEqual(row.counts['C'], 2)
                
            if row.position == 35:
                self.assertEqual(row.counts['C'], 1)
                self.assertEqual(row.counts['T'], 2)
    
    def test_base_qualities(self):
        self._check_base_qual(0, 3)
        self._check_base_qual(10, 3)
        self._check_base_qual(11, 2)
        
        self._check_base_qual(19, 2)
        self._check_base_qual(20, 2)
        self._check_base_qual(21, 1)
        
        self._check_base_qual(30, 1)
        
        self._check_base_qual(40, 0)
        
    def _check_base_qual(self, base_qual, expected_counts):
        reader = self.get_counter('base_quality.bam', min_base_qual=base_qual)        
        row = reader.iter_ref('1').next()        
        self.assertEqual(row.counts['A'], expected_counts)

    def test_map_qualities(self):
        self._check_map_qual(0, 4)
        self._check_map_qual(10, 3)
        self._check_map_qual(11, 2)
        
        self._check_map_qual(19, 2)
        self._check_map_qual(20, 2)
        self._check_map_qual(21, 1)
        
        self._check_map_qual(30, 1)
        
        self._check_map_qual(40, 0)
        
    def _check_map_qual(self, map_qual, expected_counts):
        reader = self.get_counter('map_quality.bam', min_map_qual=map_qual)        
        row = reader.iter_ref('1').next()        
        self.assertEqual(row.counts['A'], expected_counts)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
