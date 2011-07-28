'''
Suite of tests for the JointBinaryBaseCounter class.

Created on 2011-06-30

@author: Andrew Roth
'''
import unittest

import pysam

from joint_snv_mix.counters.joint_binary_quality_counter import JointBinaryQualityCounter


class Test(unittest.TestCase):
    def get_counter(self):
        normal_file_name = "data/paired/normal.bam"
        
        tumour_file_name = "data/paired/tumour.bam"
        
        reference_file_name = "data/paired/reference.fasta"
        
        ref_genome = pysam.Fastafile(reference_file_name)
        
        normal_bam = pysam.Samfile(normal_file_name, 'rb')
        tumour_bam = pysam.Samfile(tumour_file_name, 'rb')
             
        counter = JointBinaryQualityCounter(normal_bam, tumour_bam, ref_genome)

        return counter
    
    def check_position(self, row, ref, position, ref_base, non_ref_base):
        self.assertEqual(row.ref, ref)
        self.assertEqual(row.position, position)
        self.assertEqual(row.ref_base, ref_base)
        self.assertEqual(row.non_ref_base, non_ref_base)
    
    def test_str(self):
        counter = self.get_counter()
        iter = counter.iter_ref('1')
        
        row = iter.next()
        
        s = str(row)
        
        s_parts = s.split("\t")
        
        self.assertEqual(s_parts[0], '1')
        self.assertEqual(int(s_parts[1]), 1)
        self.assertEqual(s_parts[2], 'T')
        self.assertEqual(s_parts[3], 'N')
        
    def test_iter_matches_row(self):
        counter = self.get_counter()        
        iter = counter.iter_ref('1')
        
        for row in iter:
            self.assertEqual(row.position, iter.position)
            self.assertEqual(row.ref, iter.ref)
    
    def test_simple_counts(self):       
        counter = self.get_counter()
        
        iter = counter.iter_ref('1')
        
        row = iter.next()        
        self.check_position(row, '1', 1, 'T', 'N')
        
        row = iter.next()
        self.check_position(row, '1', 2, 'A', 'C')
   
    def test_tumour_del(self):
        counter = self.get_counter()
        
        for row in counter.iter_ref('1'):
            print row
            if row.position == 3:
                # Check normal
                self.assertEqual(row.counts[0][0][0], 'A')
                self.assertEqual(row.counts[0][1][0], 'A')
                self.assertEqual(row.counts[0][2][0], 'A')
                
                # Check tumour
                self.assertEqual(row.counts[1][0][0], 'A')
                self.assertEqual(row.counts[1][1][0], 'A')
    
    def test_tumour_insertion(self):
        counter = self.get_counter()
        
        for row in counter.iter_ref('1'):
            if row.position == 4:
                self.assertEqual(row.counts[1][0][0], 'C')
                self.assertEqual(row.counts[1][1][0], 'C')
                self.assertEqual(row.counts[1][2][0], 'C')
            
            if row.position == 5:
                self.assertEqual(row.counts[1][0][0], 'C')
                self.assertEqual(row.counts[1][1][0], 'C')
                self.assertEqual(row.counts[1][2][0], 'C')
                

    def test_homzygous_germline(self):
        counter = self.get_counter()
        
        for row in counter.iter_ref('1'):
            if row.position == 6:
                self.assertEqual(row.ref_base, 'C')
                self.assertEqual(row.non_ref_base, 'T')
                
    def test_hetorzygous_germline(self):
        counter = self.get_counter()
        
        for row in counter.iter_ref('1'):
            if row.position == 7:
                self.assertEqual(row.ref_base, 'T')
                self.assertEqual(row.non_ref_base, 'A')
                
    def test_loh_to_variant(self):
        counter = self.get_counter()
        
        for row in counter.iter_ref('1'):
            if row.position == 8:
                self.assertEqual(row.ref_base, 'A')
                self.assertEqual(row.non_ref_base, 'C')
                
    def test_loh_to_reference(self):
        counter = self.get_counter()
        
        for row in counter.iter_ref('1'):
            if row.position == 9:
                self.assertEqual(row.ref_base, 'A')
                self.assertEqual(row.non_ref_base, 'G')          

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
