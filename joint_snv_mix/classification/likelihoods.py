'''
Created on 2011-01-18

@author: Andrew Roth
'''
import numpy as np
from joint_snv_mix.classification.utils.log_pdf import log_beta_binomial_likelihood, log_binomial_likelihood, \
    log_multinomial_likelihood
from joint_snv_mix import constants

#=======================================================================================================================
# Independent Models
#=======================================================================================================================
def independent_beta_binomial_log_likelihood(data, parameters):
    a = data.a
    b = data.b
    
    d = a + b
    
    alpha = parameters['alpha']
    beta = parameters['beta']
    
    log_likelihoods = log_beta_binomial_likelihood(a, d, alpha, beta)

    pi = parameters['pi']
    log_likelihoods = log_likelihoods + np.log(pi)
    
    return log_likelihoods

def independent_binomial_log_likelihood(data, parameters):
    a = data.a
    b = data.b
    
    d = a + b
    
    mu = parameters['mu']

    log_likelihoods = log_binomial_likelihood(a, d, mu)

    pi = parameters['pi']
    log_pi = np.log(pi)

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods

#=======================================================================================================================
# Joint Models
#=======================================================================================================================
def joint_beta_binomial_log_likelihood(data, parameters):    
    log_likelihoods = {}
    
    for genome in constants.genomes:
        a = data.a[genome]
        b = data.b[genome]
        d = a + b
        
        alpha = parameters[genome]['alpha']
        beta = parameters[genome]['beta']
    
        log_likelihoods[genome] = log_beta_binomial_likelihood(a, d, alpha, beta)

    pi = parameters['pi']

    log_likelihoods = get_joint_log_likelihoods(log_likelihoods, pi)
    
    return log_likelihoods

def joint_binomial_log_likelihood(data, parameters):
    log_likelihoods = {}
    
    for genome in constants.genomes:
        a = data.a[genome]
        b = data.b[genome]
        d = a + b
        
        mu = parameters[genome]['mu']
    
        log_likelihoods[genome] = log_binomial_likelihood(a, d, mu)

    pi = parameters['pi']

    log_likelihoods = get_joint_log_likelihoods(log_likelihoods, pi)
    
    return log_likelihoods

def get_joint_log_likelihoods(log_likelihoods, pi):
    normal_log_likelihoods = log_likelihoods['normal']
    tumour_log_likelihoods = log_likelihoods['tumour']
    
    normal_nclass = normal_log_likelihoods.shape[1]
    column_shape = (normal_log_likelihoods[:, 0].size, 1)

    log_likelihoods = np.hstack([normal_log_likelihoods[:, i].reshape(column_shape) + tumour_log_likelihoods
                                  for i in range(normal_nclass)])
    
    log_likelihoods = log_likelihoods + np.log(pi)

    return log_likelihoods

#=======================================================================================================================
# Multinomial
#=======================================================================================================================
def joint_multinomial_log_likelihood(data, parameters):
    log_likelihoods = {}
    
    for genome in constants.genomes:
        counts = data.counts[genome]        
        rho = parameters[genome]['rho']
    
        log_likelihoods[genome] = log_multinomial_likelihood(counts, rho)

    pi = parameters['pi']
    
    log_likelihoods = get_joint_log_likelihoods(log_likelihoods, pi)
    
    return log_likelihoods
