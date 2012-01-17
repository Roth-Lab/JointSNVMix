'''
Created on 2011-06-21

@author: andrew
'''
import unittest

from joint_snv_mix.samtools import BamFile
from joint_snv_mix.counters.base_counter import BaseCounter

class Test(unittest.TestCase):
    def get_counter(self, bam_file_name, min_base_qual=10, min_map_qual=10):
        bam_file = BamFile(bam_file_name)        
        counter = BaseCounter(bam_file, min_base_qual, min_map_qual)

        return counter
    
    def test_iter_matches_row(self):
        counter = self.get_counter('data/plain.bam', 0, 0)        
        iter = counter.iter_ref('1')
        
        for row in iter:
            self.assertEqual(row.position, iter.position)
            self.assertEqual(row.ref, iter.ref)
    
    def test_str(self):
        counter = self.get_counter('data/plain.bam')
        iter = counter.iter_ref('1')
        
        row = iter.next()
        
        s = str(row)
        
        s_parts = s.split("\t")
        
        self.assertEqual(s_parts[0], '1')
        self.assertEqual(int(s_parts[1]), 33)
        
    
    def test_refs(self):
        reader = self.get_counter('data/plain.bam')
        
        self.assertTrue('1' in reader.refs)
        
    def test_iter_row(self):
        reader = self.get_counter('data/plain.bam')
        
        iter = reader.iter_ref('1')
        
        row_1 = iter.next()        
        self.assertEqual(row_1.position, 33)
        self.assertEqual(row_1.counts[0], 1)
        self.assertEqual(sum(row_1.counts), 1)
        self.assertEqual(row_1.depth, 1)
    
    def test_deletion_correct(self):
        reader = self.get_counter('data/deletion.bam')
        
        for row in reader.iter_ref('1'):
            if row.position == 34:
                self.assertEqual(sum(row.counts), 0)
                self.assertEqual(row.depth, 0)
                
            if row.position == 35:
                self.assertEqual(row.counts[1], 1)
                self.assertEqual(row.counts[2], 2)

    def test_insertion_correct(self):
        reader = self.get_counter('data/insertion.bam')
        
        for row in reader.iter_ref('1'):
            if row.position == 34:
                self.assertEqual(row.counts[1], 2)
                
            if row.position == 35:
                self.assertEqual(row.counts[1], 1)
                self.assertEqual(row.counts[3], 2)
    
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
        reader = self.get_counter('data/base_quality.bam', min_base_qual=base_qual)        
        row = reader.iter_ref('1').next()        
        self.assertEqual(row.counts[0], expected_counts)

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
        reader = self.get_counter('data/map_quality.bam', min_map_qual=map_qual)        
        row = reader.iter_ref('1').next()        
        self.assertEqual(row.counts[0], expected_counts)
    
    def test_lower_case(self):
        '''
        Make sure that lower case in bam files does not affect counts.
        '''
        counter = self.get_counter('data/lower_case.bam', 0, 0)        
        iter = counter.iter_ref('1')
        
        row = iter.next()
        
        self.assertEqual(row.counts[0], 1)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
