'''
Created on 2010-11-18

@author: Andrew Roth
'''
import numpy as np

from scipy.stats import binom, poisson


def draw_easy_sample(n):
    mean_depth = 1000
    parameters = get_easy_parameters()
    
    return draw_joint_sample(parameters, n, mean_depth)
    
def draw_joint_sample(parameters, n, mean_depth):    
    class_params = []
    
    for i in range(3):
        for j in range(3):
            class_params.append([ 
                                 parameters['normal']['mu'][i], parameters['tumour']['mu'][j] 
                                 ])
    
    multinomial_draw = np.random.multinomial(1, parameters['pi'], n)

    labels = np.argmax(multinomial_draw == 1, axis=1)
    
    a_b = []

    for label in labels:
        a_b.append(class_params[label])
    
    a_b = np.array(a_b)

    x_N = draw_binomial_sample(a_b[:, 0], mean_depth, n)
    x_T = draw_binomial_sample(a_b[:, 1], mean_depth, n)
        
    x = np.hstack((x_N, x_T))
    
    return x, labels

def draw_binomial_sample(mu, mean_depth, n):
    d = poisson.rvs(mean_depth, size=n)
    
    while np.any(d == 0):
        n = np.where(d == 0)[0].size
        
        d[d == 0] = poisson.rvs(mean_depth, size=n)
    
    a = binom.rvs(d, mu)
    
    b = d - a
    
    x = np.column_stack((a, b))
    
    return x

def get_easy_parameters():
    parameters = {}
    
    parameters['normal'] = {}
    parameters['tumour'] = {}
    
    parameters['normal']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)
    parameters['tumour']['mu'] = np.array([0.999, 0.5, 0.001], dtype=np.float)

    parameters['pi'] = np.ones((9,))
    parameters['pi'] = parameters['pi'] / parameters['pi'].sum()
    
    
    return parameters
