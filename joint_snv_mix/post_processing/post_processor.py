'''
Created on 2012-01-28

@author: Andrew Roth
'''
from joint_snv_mix.post_processing.feature_extractor import FeatureExtractor

import cpickle
import xalglib

class PostProcessor(object):
    def __init__(self, normal_bam, tumour_bam, rf_pickle_file):
        self._feature_extractor = FeatureExtractor(normal_bam, tumour_bam)
        
        self._rf = cpickle.load(open(rf_pickle_file, 'rb')) 

    def get_somatic_prob(self, row):
        features = self._feature_extractor.get_features(row)
        
        y = xalglib.dfprocessi(self._rf, features.values())
        
        return y[0]
        
        
