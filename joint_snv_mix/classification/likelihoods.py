'''
Created on 2011-01-18

@author: Andrew Roth
'''
import numpy as np
from joint_snv_mix.classification.utils.log_pdf import log_beta_binomial_likelihood, log_binomial_likelihood
from joint_snv_mix import constants

#=======================================================================================================================
# Independent Models
#=======================================================================================================================
def independent_beta_binomial_log_likelihood( data, parameters ):
    a = data.a
    b = data.b
    
    d = a + b
    
    alpha = parameters['alpha']
    beta = parameters['beta']
    
    log_likelihoods = log_beta_binomial_likelihood( a, d, alpha, beta )

    pi = parameters['pi']
    log_likelihoods = log_likelihoods + np.log( pi )
    
    return log_likelihoods

def independent_binomial_log_likelihood( data, parameters ):
    a = data.a
    b = data.b
    
    d = a + b
    
    mu = parameters['mu']

    log_likelihoods = log_binomial_likelihood( a, d, mu )

    pi = parameters['pi']
    log_pi = np.log( pi )

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods

#=======================================================================================================================
# Joint Models
#=======================================================================================================================
def joint_beta_binomial_log_likelihood( data, parameters ):    
    log_likelihoods = {}
    
    for genome in constants.genomes:
        a = data.a[genome]
        b = data.b[genome]
        d = a + b
        
        alpha = parameters[genome]['alpha']
        beta = parameters[genome]['beta']
    
        log_likelihoods[genome] = log_beta_binomial_likelihood( a, d, alpha, beta )

    pi = parameters['pi']

    log_likelihoods = get_joint_log_likelihoods( log_likelihoods, pi )
    
    return log_likelihoods

def joint_binomial_log_likelihood( data, parameters ):
    a_1 = data.a[0]
    a_2 = data.a[1]
    
    d_1 = data.a[0] + data.b[0]
    d_2 = data.a[1] + data.b[1]
    
    mu_1 = parameters['mu'][0]
    mu_2 = parameters['mu'][1]

    
    normal_log_likelihoods = log_binomial_likelihood( a_1, d_1, mu_1 )
    tumour_log_likelihoods = log_binomial_likelihood( a_2, d_2, mu_2 )

    pi = parameters['pi']

    log_likelihoods = get_joint_log_likelihoods( normal_log_likelihoods, tumour_log_likelihoods, pi )
    
    return log_likelihoods

def get_joint_log_likelihoods( log_likelihoods, pi ):
    normal_log_likelihoods = log_likelihoods['normal']
    tumour_log_likelihoods = log_likelihoods['tumour']
    
    column_shape = ( normal_log_likelihoods[:, 0].size, 1 )

    log_likelihoods = np.hstack( ( 
                                 normal_log_likelihoods[:, 0].reshape( column_shape ) + tumour_log_likelihoods ,
                                 normal_log_likelihoods[:, 1].reshape( column_shape ) + tumour_log_likelihoods ,
                                 normal_log_likelihoods[:, 2].reshape( column_shape ) + tumour_log_likelihoods
                                 ) )
    
    log_likelihoods = log_likelihoods + np.log( pi )

    return log_likelihoods
