'''
Created on 2011-03-31

@author: andrew
'''
import unittest

import numpy as np

from joint_snv_mix import constants
from joint_snv_mix.classification.utils.data import JointData
from joint_snv_mix.classification.fisher_test import IndependentFisherModel

class Test(unittest.TestCase):
    def setUp(self):
        self.counts = np.array([
                               [   1000, 0, 1000, 0],
                               [   1000, 0, 500, 500],
                               [   1000, 0, 0, 1000],
                               [   500, 500, 1000, 0],
                               [   500, 500, 500, 500],
                               [   500, 500, 0, 1000],
                               [   0, 1000, 1000, 0],
                               [   0, 1000, 500, 500],
                               [   0, 1000, 0, 1000],
                              ])
        
        self.data = JointData(self.counts)
        
        self.labels = np.arange(9)
        
        self.model = IndependentFisherModel()
    
    def test_easy_classification(self):
        '''
        Test simple calls that all methods should make correctly.
        '''
        predicted_labels = self.model.classify(self.data)
                
        for pred_label, label in zip(predicted_labels, self.labels):
            err_str = "Genotype {0} is predicted as {1}.".format(constants.joint_genotypes[label],
                                                                 constants.joint_genotypes[pred_label])
            
            self.assertEqual(pred_label, label, err_str)
            
    def test_noisy_classification(self):
        '''
        Test simple calls with a small amount of noise added.
        '''
        min_noise = 0
        max_noise = 30
        
        noise = np.random.randint(min_noise, max_noise, size=(9, 4))
        
        counts = self.counts + noise
        data = JointData(counts)
        
        predicted_labels = self.model.classify(data)
                
        for pred_label, label in zip(predicted_labels, self.labels):
            err_str = "Genotype {0} is predicted as {1}.".format(constants.joint_genotypes[label],
                                                                 constants.joint_genotypes[pred_label])
            
            self.assertEqual(pred_label, label, err_str)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
