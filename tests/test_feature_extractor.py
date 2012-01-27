'''
Created on 2012-01-26

@author: Andrew Roth
'''
import unittest

from joint_snv_mix.counter import JointBinaryCounter
from joint_snv_mix.samtools import BamFile, FastaFile

from joint_snv_mix.post_processing.feature_extractor import FeatureExtractor

class Test(unittest.TestCase):
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
    
    def test_feature_extractor(self):
        counter = self.get_counter(0, 0)
        
        iter = counter.get_ref_iterator('1')
        
        extractor = FeatureExtractor(self._normal_bam, self._tumour_bam)
        
        row = iter.next()
        print extractor.get_features(row)
        
        row = iter.next()
        print extractor.get_features(row)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
