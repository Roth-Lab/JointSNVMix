'''
Created on 2012-01-17

@author: andrew
'''
import unittest

from joint_snv_mix.samtools.bam import BamFile


class Test(unittest.TestCase):
    def test_references(self):
        bam = BamFile('data/multiple_chromosome.bam')
        
        self.assertSetEqual(set(('1', '2')), set(bam.references))
        
    def test_references(self):
        bam = BamFile('data/deletion.bam')
        
        iter = bam.get_counts_iterator('1')
        
        for x in iter:
            print x  


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
