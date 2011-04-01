'''
Created on 2011-03-31

@author: andrew
'''
import unittest

import numpy as np

from simulation.binomial import draw_easy_sample
from joint_snv_mix.classification.utils.data import JointData
from joint_snv_mix.classification.fisher import IndependentFisherModel

class Test(unittest.TestCase):
    def test_classification_error_is_sane(self):
        sample, labels = draw_easy_sample()
        
        marg_labels = self.convert_joint_genotypes_to_marginal_genotypes(labels)
        
        data = JointData(sample)
        
        model = IndependentFisherModel()

        predicted_labels = model.classify(data)
        
        
        self.assertEqual(predicted_labels == 2, marg_labels == 2)
        
    def convert_joint_genotypes_to_marginal_genotypes(self, labels):
        marginal_labels = np.zeros(labels.shape)
        
        ref_indices = (labels == 0)
        som_indices = np.logical_or(labels == 1, labels == 2)
        ger_indices = np.logical_or(labels == 4, labels == 8)
        loh_indices = np.logical_or(labels == 3, labels == 5)
        unk_indices = np.logical_or(labels == 6, labels == 7)
        
        marginal_labels[ref_indices] = 0
        marginal_labels[ger_indices] = 1
        marginal_labels[som_indices] = 2
        marginal_labels[loh_indices] = 3
        marginal_labels[unk_indices] = 4
        
        return marginal_labels


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
