'''
Created on 2011-09-07

@author: andrew
'''
import unittest

from pysam import Samfile

from joint_snv_mix.counters.base_counter import BaseCounter
from joint_snv_mix.counters.positions_counter import PositionsCounter

class Test(unittest.TestCase):
    def get_positions_counter(self, positions_file):
        bam_file = Samfile("data/positions/plain.bam", 'rb')
        counter = BaseCounter(bam_file, 0, 0)
        
        pos_counter = PositionsCounter(counter, positions_file)
        
        return pos_counter

    def test_iterates_over_simple_positions(self):
        positions_file = "data/positions/simple.txt"
        
        counter = self.get_positions_counter(positions_file)
        
        iter = counter.iter_ref('1')
        
        row = iter.next()       
        self.assertEqual(row.position, 33)
        
        row = iter.next()       
        self.assertEqual(row.position, 38)
        
    def test_iterates_positions_with_extra_chrom(self):
        positions_file = "data/positions/missing_chrom.txt"
        
        counter = self.get_positions_counter(positions_file)
        
        iter = counter.iter_ref('3')
        
        row = iter.next()
        self.assertEqual(row.position, 100)
        
        row = iter.next()       
        self.assertEqual(row.position, 102)        
    
    def test_raises_exception_when_ref_not_present(self):
        positions_file = "data/positions/simple.txt"
        
        counter = self.get_positions_counter(positions_file)
        
        self.assertRaises(Exception, counter.iter_ref, "X")
        
    def test_iterates_over_pos_file_with_extra_chrom(self):
        '''
        Regression test to make sure it doesn't error out when positions file has chromosomes not present in bam file.
        '''
        positions_file = "data/positions/missing_chrom.txt"
        
        counter = self.get_positions_counter(positions_file)
        
        self.assertEqual(('1', '3'), counter.refs)
        
        for ref in counter.refs:
            for row in counter.iter_ref(ref):
                pass

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
