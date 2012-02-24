'''
Created on 2012-01-28

@author: Andrew Roth
'''
from joint_snv_mix.post_processing.feature_extractor import FeatureExtractor

import bz2
import xalglib

class PostProcessor(object):
    def __init__(self, normal_bam, tumour_bam, rf_pickle_file):
        self._feature_extractor = FeatureExtractor(normal_bam, tumour_bam)
        
        import os
        print os.getcwd()
        
        self._rf = xalglib.dfunserialize(bz2.BZ2File(rf_pickle_file, 'r').read()) 

    def get_somatic_prob(self, row):
        features = self._feature_extractor.get_features(row)
        
        y = xalglib.dfprocessi(self._rf, features.values())
        
        return y[1]
        
        
