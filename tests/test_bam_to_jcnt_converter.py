'''
Created on 2011-03-10

@author: andrew
'''
import unittest

import pysam

from joint_snv_mix.pre_processing.bam_to_jcnt import BamToJcntConverter, JointPileupIterator


class Test(unittest.TestCase):
    def setUp(self):
        self.normal_bam = pysam.Samfile('data/normal.bam', 'rb')
        self.tumour_bam = pysam.Samfile('data/tumour.bam', 'rb')
        
        self.ref_genome_fasta = FastaFileFake()        
    
    def test_get_counts_no_non_ref(self):
        converter = BamToJcntConverter(self.normal_bam, self.tumour_bam, self.ref_genome_fasta)
        
        ref_base = 'A'
        bases = ['A', 'A', 'A', 'C', 'C', 'T', 'A']
        
        non_ref_base, counts = converter._get_counts(ref_base, bases)
        
        self.assertEqual(4, counts[0])
        self.assertEqual('C', non_ref_base)
        self.assertEqual(2, counts[1])

    def test_get_counts_with_non_ref(self):
        converter = BamToJcntConverter(self.normal_bam, self.tumour_bam, self.ref_genome_fasta)
        
        ref_base = 'A'
        bases = ['A', 'A', 'A', 'C', 'C', 'T', 'A']
        
        non_ref_base, counts = converter._get_counts(ref_base, bases, non_ref_base="T")
        
        self.assertEqual(4, counts[0])
        self.assertEqual('T', non_ref_base)
        self.assertEqual(1, counts[1])
        
    def test_min_depth(self):
        converter = BamToJcntConverter(self.normal_bam,
                                       self.tumour_bam,
                                       self.ref_genome_fasta,
                                       min_depth=10,
                                       min_bqual=0,
                                       min_mqual=0)
        
        ref = 'Y'
        start = 2905239
        stop = 2905239
        
        iter = JointPileupIterator(
                                   self.normal_bam.pileup(ref, start, stop),
                                   self.tumour_bam.pileup(ref, start, stop)
                                   )
        
        for normal_col, tumour_col in iter:
            pos = normal_col.pos
            
            if pos == 2905239:            
                entry = converter._get_jcnt_entry(normal_col, tumour_col, pos, 'A')
            
                self.assertEqual(entry, None)
                
    def test_min_bqual(self):
        converter = BamToJcntConverter(self.normal_bam,
                                       self.tumour_bam,
                                       self.ref_genome_fasta,
                                       min_depth=0,
                                       min_bqual=21,
                                       min_mqual=0)
        
        ref = 'Y'
        start = 2905239
        stop = 2905239
        
        iter = JointPileupIterator(
                                   self.normal_bam.pileup(ref, start, stop),
                                   self.tumour_bam.pileup(ref, start, stop)
                                   )
        
        for normal_col, tumour_col in iter:
            pos = normal_col.pos
            
            if pos == 2905239:            
                entry = converter._get_jcnt_entry(normal_col, tumour_col, pos, 'A')
            
                self.assertEqual(entry[4], 7)
                self.assertEqual(entry[5], 0)
                self.assertEqual(entry[6], 9)
                self.assertEqual(entry[7], 0)
                
    def test_min_mqual(self):
        converter = BamToJcntConverter(self.normal_bam,
                                       self.tumour_bam,
                                       self.ref_genome_fasta,
                                       min_depth=0,
                                       min_bqual=0,
                                       min_mqual=27)
        
        ref = 'Y'
        start = 2905239
        stop = 2905239
        
        iter = JointPileupIterator(
                                   self.normal_bam.pileup(ref, start, stop),
                                   self.tumour_bam.pileup(ref, start, stop)
                                   )
        
        for normal_col, tumour_col in iter:
            pos = normal_col.pos
            
            if pos == 2905239:            
                entry = converter._get_jcnt_entry(normal_col, tumour_col, pos, 'A')
            
                self.assertEqual(entry[4], 3)
                self.assertEqual(entry[5], 0)
                self.assertEqual(entry[6], 2)
                self.assertEqual(entry[7], 0)        

class BamFileFake:
    def __init__(self):
        self.references = ['1', '2'] 

class FastaFileFake:
    def __init__(self):
        pass

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
