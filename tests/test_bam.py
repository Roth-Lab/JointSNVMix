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


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
