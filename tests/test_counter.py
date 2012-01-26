'''
Suite of tests for the JointBinaryBaseCounter class.

Created on 2011-06-30

@author: Andrew Roth
'''
from __future__ import division

import unittest

from joint_snv_mix.counter import JointBinaryCounter

from joint_snv_mix.samtools.bam import BamFile
from joint_snv_mix.samtools.fasta import FastaFile

class CounterTester(object):   
    def check_position(self, row, ref, position, ref_base, var_base):
        self.assertEqual(row.ref, ref)
        self.assertEqual(row.position, position)
        self.assertEqual(row.ref_base, ref_base)
        self.assertEqual(row.var_base, var_base)
    
    def check_counts(self, row, counts):
        self.assertTupleEqual(row.counts, counts)

    def test_str(self):
        counter = self.get_counter(0, 0)
        iter = counter.get_ref_iterator('1')
        
        row = iter.next()
        
        s = str(row)

        s_parts = s.split("\t")
        
        self.assertEqual(s_parts[0], '1')
        self.assertEqual(int(s_parts[1]), 1)
        self.assertEqual(s_parts[2], 'T')
        self.assertEqual(s_parts[3], 'N')
    
    def test_simple_counts(self):       
        counter = self.get_counter(0, 0)
        
        iter = counter.get_ref_iterator('1')
        
        row = iter.next()        
        self.check_position(row, '1', 1, 'T', 'N')
        
        row = iter.next()
        self.check_position(row, '1', 2, 'A', 'C')
    
    def test_tumour_del(self):
        counter = self.get_counter(0, 0)
        
        for row in counter.get_ref_iterator('1'):
            if row.position == 3:
                self.check_counts(row, (3, 0, 2, 0))
    
    def test_tumour_insertion(self):
        counter = self.get_counter(0, 0)
        
        for row in counter.get_ref_iterator('1'):
            if row.position in [4, 5]:
                self.check_counts(row, (3, 0, 3, 0))

    def test_homzygous_germline(self):
        counter = self.get_counter(0, 0)
        
        for row in counter.get_ref_iterator('1'):
            if row.position == 6:
                self.assertEqual(row.ref_base, 'C')
                self.assertEqual(row.var_base, 'T')
                
                self.check_counts(row, (0, 3, 0, 3))
                
    def test_hetorzygous_germline(self):
        counter = self.get_counter(0, 0)
        
        for row in counter.get_ref_iterator('1'):
            if row.position == 7:
                self.assertEqual(row.ref_base, 'T')
                self.assertEqual(row.var_base, 'A')
                
                self.check_counts(row, (1, 2, 2, 1))
                
    def test_loh_to_variant(self):
        counter = self.get_counter(0, 0)
        
        for row in counter.get_ref_iterator('1'):
            if row.position == 8:
                self.assertEqual(row.ref_base, 'A')
                self.assertEqual(row.var_base, 'C')
                
                self.check_counts(row, (1, 2, 0, 3))
                
    def test_loh_to_reference(self):
        counter = self.get_counter(0, 0)
        
        for row in counter.get_ref_iterator('1'):
            if row.position == 9:
                self.assertEqual(row.ref_base, 'A')
                self.assertEqual(row.var_base, 'G')
                
                self.check_counts(row, (2, 1, 3, 0))                

    def test_lower_case(self):
        '''
        Make sure lower case letters in BAM file cause no issues.
        '''
        counter = self.get_counter(0, 0)
        
        for row in counter.get_ref_iterator('1'):
            if row.position == 10:
                self.assertEqual(row.ref_base, 'C')
                self.assertEqual(row.var_base, 'G')
                
                self.check_counts(row, (3, 0, 2, 1))

class TestBaseCounter(CounterTester, unittest.TestCase):
    def setUp(self):
        normal_file_name = "data/paired/normal.bam"
        
        tumour_file_name = "data/paired/tumour.bam"
        
        reference_file_name = "data/paired/reference.fasta"
        
        self._ref_genome = FastaFile(reference_file_name)
        
        self._normal_bam = BamFile(normal_file_name)
        self._tumour_bam = BamFile(tumour_file_name)
            
    def get_counter(self, min_base_qual=10, min_map_qual=10):
        counter = JointBinaryCounter(self._normal_bam, self._tumour_bam, self._ref_genome, min_base_qual, min_map_qual)

        return counter
    
    def test_map_qualities(self):
        counter = self.get_counter(0, 0)        
        iter = counter.get_ref_iterator('1')
        row = iter.next()
        
        self.check_counts(row, (3, 0, 3, 0))
        
        counter = self.get_counter(0, 10)        
        iter = counter.get_ref_iterator('1')
        row = iter.next()        
        self.check_counts(row, (3, 0, 3, 0))
        
        counter = self.get_counter(0, 20)        
        iter = counter.get_ref_iterator('1')
        row = iter.next()
        self.check_counts(row, (2, 0, 2, 0))
        
        counter = self.get_counter(0, 30)        
        iter = counter.get_ref_iterator('1')
        row = iter.next()
        self.check_counts(row, (1, 0, 1, 0))    

class TestQualitiesCounter(CounterTester, unittest.TestCase):
    def setUp(self):
        normal_file_name = "data/paired/normal.bam"
        
        tumour_file_name = "data/paired/tumour.bam"
        
        reference_file_name = "data/paired/reference.fasta"
        
        self._ref_genome = FastaFile(reference_file_name)
        
        self._normal_bam = BamFile(normal_file_name)
        self._tumour_bam = BamFile(tumour_file_name)
            
    def get_counter(self, min_base_qual=10, min_map_qual=10):
        counter = JointBinaryCounter(self._normal_bam, self._tumour_bam, self._ref_genome, min_base_qual, min_map_qual,
                                     qualities=1)

        return counter
    
    def test_quality_probabilities(self):
        counter = self.get_counter(0, 0)        
        iter = counter.get_ref_iterator('1')
        row = iter.next()
        
        expected_map_probs = [self.convert_phred_to_probabilities(x) for x in [10, 20, 30]]        
        self.assertSequenceEqual(row.data.normal_mapping_qualities, expected_map_probs)
        
        expected_map_probs = [self.convert_phred_to_probabilities(x) for x in [10, 25, 30]]        
        self.assertSequenceEqual(row.data.tumour_mapping_qualities, expected_map_probs)
        
        row = iter.next()
               
        expected_base_probs = [self.convert_phred_to_probabilities(x) for x in [27, 27, 27]]        
        self.assertSequenceEqual(row.data.normal_base_qualities, expected_base_probs)
        
        expected_base_probs = [self.convert_phred_to_probabilities(x) for x in [27, 27, 27]]
        expected_base_probs[2] = (1 - expected_base_probs[2]) / 3
        
        for p, expected_p in zip(row.data.tumour_base_qualities, expected_base_probs):    
            self.assertAlmostEqual(p, expected_p, 5)
                                                            
    def convert_phred_to_probabilities(self, phred_qual):
        return 1 - 10 ** (-phred_qual / 10)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
