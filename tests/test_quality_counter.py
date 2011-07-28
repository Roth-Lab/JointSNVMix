'''
Created on 2011-06-21

@author: andrew
'''
import unittest

import pysam

from joint_snv_mix.counters.quality_counter import QualityCounter

class Test(unittest.TestCase):
    def get_counter(self, bam_file_name):
        bam_file = pysam.Samfile(bam_file_name, 'rb')        
        counter = QualityCounter(bam_file)

        return counter
    
    def test_iter_matches_row(self):
        counter = self.get_counter('data/plain.bam')        
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
        
        iter.next()
        iter.next()
        row_1 = iter.next()     
        self.assertEqual(row_1.position, 35)
        expect_row = (
                      'C',
                      self._get_prob_from_phred(ord('<') - 33),
                      self._get_prob_from_phred(20)
                      )
        self.assertEqual(row_1.counts[0], expect_row)
        self.assertEqual(row_1.depth, 3)

    def _get_prob_from_phred(self, qual):
        return 1. - pow(10, -qual / 10.) 
    
    def test_deletion_correct(self):
        reader = self.get_counter('data/deletion.bam')
        
        for row in reader.iter_ref('1'):
            if row.position == 34:
                self.assertEqual(row.depth, 0)
                
            if row.position == 35:
                self.assertEqual(row.counts[0][0], ('G'))
                self.assertEqual(row.counts[1][0], ('G'))
                self.assertEqual(row.counts[2][0], ('C'))

    def test_insertion_correct(self):
        reader = self.get_counter('data/insertion.bam')
        
        for row in reader.iter_ref('1'):
            if row.position == 34:
                self.assertEqual(row.counts[0][0], ('C'))
                self.assertEqual(row.counts[1][0], ('C'))
                                
            if row.position == 35:
                self.assertEqual(row.counts[0][0], ('T'))
                self.assertEqual(row.counts[1][0], ('T'))
                self.assertEqual(row.counts[2][0], ('C'))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
