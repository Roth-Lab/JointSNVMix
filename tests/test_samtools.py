'''
Created on 2012-01-17

@author: andrew
'''
import unittest

from joint_snv_mix.samtools.bam import BamFile
from joint_snv_mix.samtools.fasta import FastaFile


class Test(unittest.TestCase):
    def test_references(self):
        bam = BamFile('data/pileup.bam')
        
        self.assertSetEqual(set(('1', '2', 'chrX')), set(bam.references))
    
    def test_invalid_reference(self):
        bam = BamFile('data/pileup.bam')
        
        self.assertRaises(Exception, bam.get_pileup_iterator, 'X')
        
    def test_pileup(self):
        bam = BamFile('data/pileup.bam')
        
        pileup_iter = bam.get_pileup_iterator('1')
        
        pileup_column = pileup_iter.next()
        self.assertEqual(pileup_column.position, 33)      
        self.assertListEqual(pileup_column.base_quals, [0])
        self.assertListEqual(pileup_column.map_quals, [10])
        
        pileup_column = pileup_iter.next()
        self.assertEqual(pileup_column.position, 34)        
        self.assertListEqual(pileup_column.base_quals, [0, 30])
        self.assertListEqual(pileup_column.map_quals, [10, 20])

        pileup_column = pileup_iter.next()
        self.assertEqual(pileup_column.position, 35)      
        self.assertListEqual(pileup_column.base_quals, [0, 30, 93])
        self.assertListEqual(pileup_column.map_quals, [10, 20, 100]) 
    
    def test_fasta(self):
        fasta = FastaFile('data/paired/reference.fasta')
        
        self.assertEqual(fasta.get_base('1', 1), 'T')
        self.assertEqual(fasta.get_base('1', 2), 'A')
        self.assertEqual(fasta.get_base('1', 3), 'A')
        self.assertEqual(fasta.get_base('1', 4), 'C')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
