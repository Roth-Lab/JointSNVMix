'''
Created on 2011-04-01

@author: andrew
'''
import unittest

import numpy as np
import numpy.testing as npt

from scipy.stats import binom, poisson

from simulation.binomial import get_easy_parameters, draw_easy_sample
from joint_snv_mix.classification.utils.data import JointData
from joint_snv_mix.classification.joint_binomial import JointBinomialModel


class JointBinomialTrainingTest(unittest.TestCase):
    def setUp(self):
        self.model = JointBinomialModel()
        
    def test_ml_training(self):
        '''
        Check whether maximum likelihood estimation is within some threshold of the true values sampled from the model.
        '''
        real_parameters = self._get_parameters()
        
        self.check_parameter_estimation(real_parameters)
        
        real_parameters = self._get_skew_parameters()
        
        self.check_parameter_estimation(real_parameters)

    
    def check_parameter_estimation(self, real_parameters):
        decimal_agreement = 3
        
        sample_size = int(1e6)
        
        mean_depth = 1000
        
        simulator = JointBinomialSimulator(real_parameters, mean_depth)
        
        sample, labels = simulator.draw_joint_sample(sample_size)
        
        data = JointData(sample)
        priors = self._get_priors()
        pred_parameters = self.model.train(data, priors, 1000, 1e-6)
        
        npt.assert_array_almost_equal(real_parameters['normal']['mu'],
                                      pred_parameters['normal']['mu'],
                                      decimal_agreement)
        
        npt.assert_array_almost_equal(real_parameters['tumour']['mu'],
                                      pred_parameters['tumour']['mu'],
                                      decimal_agreement)

        npt.assert_array_almost_equal(real_parameters['pi'],
                                      pred_parameters['pi'],
                                      decimal_agreement)    
    
    def _get_priors(self):
        priors = {}
        
        priors['kappa'] = np.array([1.] * 9)
        
        priors['normal'] = {}
        priors['tumour'] = {}
        
        priors['normal']['mu'] = {}
        priors['tumour']['mu'] = {}
        
        priors['normal']['mu']['alpha'] = np.array([1, 1, 1])
        priors['tumour']['mu']['alpha'] = np.array([1, 1, 1])
        
        priors['normal']['mu']['beta'] = np.array([1, 1, 1])
        priors['tumour']['mu']['beta'] = np.array([1, 1, 1])
        
        return priors
    
    def _get_parameters(self):
        parameters = {}
        
        parameters['normal'] = {}
        parameters['tumour'] = {}
        
        parameters['normal']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
        parameters['tumour']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
    
        parameters['pi'] = np.ones((9,))
        parameters['pi'] = parameters['pi'] / parameters['pi'].sum()
        
        
        return parameters
    
    def _get_skew_parameters(self):
        parameters = {}
        
        parameters['normal'] = {}
        parameters['tumour'] = {}
        
        parameters['normal']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
        parameters['tumour']['mu'] = np.array([0.999, 0.8, 0.001], dtype=np.float)
    
        parameters['pi'] = np.array((10000, 100, 100, 100, 1000, 100, 10, 10, 1000))
        parameters['pi'] = parameters['pi'] / parameters['pi'].sum()
        
        return parameters

class JointBinomialSimulator(object):
    def __init__(self, parameters, mean_depth):
        self.parameters = parameters
        self.mean_depth = mean_depth 
           
    def draw_joint_sample(self, sample_size):
        parameters = self.parameters
            
        class_params = []
        
        for i in range(3):
            for j in range(3):
                class_params.append([ 
                                     parameters['normal']['mu'][i],
                                     parameters['tumour']['mu'][j] 
                                     ])
        
        multinomial_draw = np.random.multinomial(1, parameters['pi'], sample_size)
    
        labels = np.argmax(multinomial_draw == 1, axis=1)
        
        a_b = []
    
        for label in labels:
            a_b.append(class_params[label])
        
        a_b = np.array(a_b)
    
        x_N = self.draw_binomial_sample(a_b[:, 0], sample_size)
        x_T = self.draw_binomial_sample(a_b[:, 1], sample_size)
            
        x = np.hstack((x_N, x_T))
        
        return x, labels
    
    def draw_binomial_sample(self, mu, sample_size):
        d = poisson.rvs(self.mean_depth, size=sample_size)
        
        # Keep sampling sites with depth 0 until they have counts.
        while np.any(d == 0):
            n = np.where(d == 0)[0].size
            
            d[d == 0] = poisson.rvs(self.mean_depth, size=n)
        
        a = binom.rvs(d, mu)
        
        b = d - a
        
        x = np.column_stack((a, b))
        
        return x
    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
