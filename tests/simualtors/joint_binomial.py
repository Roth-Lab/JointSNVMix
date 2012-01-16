'''
Created on 2011-08-18

@author: Andrew Roth
'''
from __future__ import division

import numpy as np

import scipy.stats as stats

class JointSnvMixSimulator(object):
    '''
    Class for drawing synthetic data from the JointSnvMixModel.
    '''
    def __init__(self, mu_N=None, mu_T=None, pi=None, normal_depth=10, tumour_depth=10):
        params = JointSnvMixParameters(mu_N, mu_T, pi)
        
        self.params = params
        self.normal_depth = normal_depth
        self.tumour_depth = tumour_depth
        
    def draw_sample(self, sample_size):
        labels = self._sample_class_labels(sample_size)
        
        mu_N = self._get_normal_density_params(labels)
        mu_T = self._get_tumour_density_params(labels)
        
        d_N = self._sample_depths(self.normal_depth, sample_size)
        d_T = self._sample_depths(self.tumour_depth, sample_size)
        
        x_N = self._draw_genome_sample(mu_N, d_N, sample_size)
        x_T = self._draw_genome_sample(mu_T, d_T, sample_size)
            
        x = np.hstack((x_N, x_T))
        
        return x, labels
    
    def _sample_class_labels(self, sample_size):
        pi = self.params.pi.flatten()
        
        multinomial_draw = np.random.multinomial(1, pi, sample_size)
        
        indices = multinomial_draw == 1
        
        labels = np.zeros((sample_size, 9))
        
        labels[indices] = 1
        
        return labels
    
    def _get_normal_density_params(self, labels):
        indices = self._get_marginal_indices(labels, "normal")
        
        mu = self.params.mu_N[:, indices].flatten()
        
        return mu
    
    def _get_tumour_density_params(self, labels):        
        indices = self._get_marginal_indices(labels, "tumour")
        
        mu = self.params.mu_T[:, indices].flatten()
        
        return mu
    
    def _get_marginal_indices(self, labels, sample):
        n = labels.shape[0]
        shape = (n, 3, 3)
        
        labels = labels.reshape(shape)
        
        if sample == "normal":
            labels = labels.sum(axis=2)
        elif sample == "tumour":        
            labels = labels.sum(axis=1)
        
        indices = np.argmax(labels, axis=1)
        
        return indices
    
    def _sample_depths(self, mean_depth, sample_size):
        # Do this to avoid zero depth and still keep mean depth.
        depths = stats.poisson.rvs(mean_depth, size=sample_size)
        
        while np.any(depths == 0):
            sample_size = np.where(depths == 0)[0].size
            
            depths[depths == 0] = stats.poisson.rvs(mean_depth, size=sample_size)
        
        return depths
    
    def _draw_genome_sample(self, mu, depth, sample_size):        
        a = stats.binom.rvs(depth, mu)        
        b = depth - a
        
        x = np.column_stack((a, b))
        
        return x

class JointSnvMixParameters(object):
    def __init__(self, mu_N=None, mu_T=None, pi=None):
        density_params_shape = (1, 3)
        class_params_shape = (1, 9)
        
        if mu_N is None:
            self.mu_N = np.array([0.99, 0.5, 0.01]).reshape(density_params_shape)
        else:
            self.mu_N = np.array(mu_N)
            
        if mu_T is None:
            self.mu_T = np.array([0.99, 0.5, 0.01]).reshape(density_params_shape)
        else:
            self.mu_T = np.array(mu_T)
            
        if pi is None:
            self.pi = np.array([1e6, 1e2, 1e2, 1e2, 1e4, 1e2, 2, 2, 1e4]).reshape(class_params_shape)
            self.pi = self.pi / self.pi.sum()
        else:
            self.pi = np.array(pi)    
    
if __name__ == "__main__":
    s = JointSnvMixSimulator()
    
    sample, labels = s.draw_sample(100000)
    
    print sample
    print labels
