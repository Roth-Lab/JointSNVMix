'''
Basic sanity tests to ensure that the classifiers can make the very easy calls.

Created on 2011-03-31

@author: Andrew Roth
'''
import unittest

import numpy as np

from joint_snv_mix import constants
from joint_snv_mix.classification.utils.data import JointData
from joint_snv_mix.classification.deterministic import IndependentFisherModel, JointFisherModel, ThresholdModel
from joint_snv_mix.classification.indep_binomial import IndependentBinomialModel
from joint_snv_mix.classification.indep import PairedIndependentModel
from joint_snv_mix.classification.joint_binomial import JointBinomialModel

class ClassifierUnitTest(object):
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
    
    def test_easy_classification(self):
        '''
        Test simple calls that all methods should make correctly.
        '''
        predicted_labels = self.classify(self.data)
        
        self.check_labels(predicted_labels)
        
    def test_noisy_classification(self):
        '''
        Test simple calls with a small amount of noise added.
        '''
        min_noise = 0
        max_noise = 50
        
        noise = np.random.randint(min_noise, max_noise, size=self.counts.shape)
        
        counts = self.counts + noise
        data = JointData(counts)
        
        predicted_labels = self.classify(data)
        
        self.check_labels(predicted_labels)
        
    def classify(self, data):
        raise NotImplemented
    
    def check_labels(self, predicted_labels):
        for pred_label, label in zip(predicted_labels, self.labels):
            err_str = "Genotype {0} is predicted as {1}.".format(constants.joint_genotypes[label],
                                                                 constants.joint_genotypes[pred_label])
            
            self.assertEqual(pred_label, label, err_str)

#===========================================================================================================
# Determinisitc Methods
#===========================================================================================================
class DeterministicClassifierUnitTest(ClassifierUnitTest):
    def classify(self, data):
        return self.model.classify(data)

class IndependentFisherUnitTest(DeterministicClassifierUnitTest, unittest.TestCase):
    def setUp(self):
        self.model = IndependentFisherModel()
        
        ClassifierUnitTest.setUp(self)

class JointFisherUnitTest(DeterministicClassifierUnitTest, unittest.TestCase):
    def setUp(self):
        self.model = JointFisherModel()
        
        ClassifierUnitTest.setUp(self)
        
class ThresholdUnitTest(DeterministicClassifierUnitTest, unittest.TestCase):
    def setUp(self):
        '''
        The threshold method only classifies reference, somatic and germline. It will fail for LOH and unknown sites.
        '''
        self.model = ThresholdModel()
        
        self.counts = np.array([
                               [   1000, 0, 1000, 0],
                               [   1000, 0, 500, 500],
                               [   1000, 0, 0, 1000],
                               [   500, 500, 500, 500],
                               [   0, 1000, 0, 1000],
                              ])
        
        self.data = JointData(self.counts)
        
        self.labels = np.array([0, 1, 2, 4, 8])

#===============================================================================================================
# Probabilistic Methods
#===============================================================================================================
class ProbabilisticClassifierUnitTest(ClassifierUnitTest):
    def classify(self, data):
        resp = self.model.classify(data, self.parameters)
        
        return np.argmax(resp, axis=1)
        
class IndependentBinomialUnitTest(ProbabilisticClassifierUnitTest, unittest.TestCase):
    def setUp(self):
        self.model = PairedIndependentModel(IndependentBinomialModel())
        
        self.parameters = self._get_parameters()
        
        ClassifierUnitTest.setUp(self)
        
    def _get_parameters(self):
        parameters = {}
        
        parameters['normal'] = {}
        parameters['tumour'] = {}
        
        parameters['normal']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
        parameters['tumour']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
        
        parameters['normal']['pi'] = np.array((1000, 100, 100), dtype=np.float)
        parameters['normal']['pi'] = parameters['normal']['pi'] / parameters['normal']['pi'].sum()
        
        parameters['tumour']['pi'] = np.array((1000, 100, 100), dtype=np.float)
        parameters['tumour']['pi'] = parameters['tumour']['pi'] / parameters['tumour']['pi'].sum()
        
        return parameters

class JointBinomialUnitTest(ProbabilisticClassifierUnitTest, unittest.TestCase):
    def setUp(self):
        self.model = JointBinomialModel()
        
        self.parameters = self._get_parameters()
        
        ClassifierUnitTest.setUp(self)
        
    def _get_parameters(self):
        parameters = {}
        
        parameters['normal'] = {}
        parameters['tumour'] = {}
        
        parameters['normal']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
        parameters['tumour']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
        
        parameters['pi'] = np.array((100000, 100, 100, 100, 1000, 100, 1, 1, 1000), dtype=np.float)
        parameters['pi'] = parameters['pi'] / parameters['pi'].sum()
        
        return parameters
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
