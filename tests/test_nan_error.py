'''
Created on 2012-01-20

@author: innovation
'''
import unittest

from joint_snv_mix.counter import JointBinaryCounter
from joint_snv_mix.models.joint_snv_mix import JointSnvMixPriors, JointSnvMixParameters, JointSnvMixModel
from joint_snv_mix.samtools.bam import BamFile
from joint_snv_mix.samtools.fasta import FastaFile
from joint_snv_mix.run import FileWriter


class Test(unittest.TestCase):
    def setUp(self):
        self.normal_bam_file = 'data/nan_error/normal.bam'
        self.tumour_bam_file = 'data/nan_error/tumour.bam'
        
        self.ref_genome_file = 'data/nan_error/ref.fasta'
        
    def get_counter(self):
        normal_bam = BamFile(self.normal_bam_file)
        tumour_bam = BamFile(self.tumour_bam_file)
        
        ref_genome = FastaFile(self.ref_genome_file)
        
        counter = JointBinaryCounter(normal_bam,
                                     tumour_bam,
                                     ref_genome,
                                     min_base_qual=0,
                                     min_map_qual=0,
                                     qualities=1)
        
        return counter
    
    def test_nan_error(self):
        counter = self.get_counter()
        
        priors = JointSnvMixPriors()
        params = JointSnvMixParameters()
        
        model = JointSnvMixModel(priors, params, model='jsm2')
        
        writer = FileWriter()
        
        for row in counter.get_ref_iterator('1'):
            if row.position == 2508:
                probs = model.predict(row.data)
                
                print probs
                
                writer.write_position(row, probs)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
