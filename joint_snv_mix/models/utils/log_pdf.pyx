'''
Created on 2011-08-04

@author: Andrew Roth
'''
#=======================================================================================================================
# SNVMix1 Code
#=======================================================================================================================
cdef double multinomial_log_likelihood(int * x,
                                       double * mu,
                                       int num_classes):
    '''
    Return the log multinomial likelihood.
    '''
    cdef int i
    cdef double log_likelihood
    
    log_likelihood = 0
    
    for i in range(num_classes):
        log_likelihood += x[i] * log(mu[i]) 

    return log_likelihood

cdef double dirichlet_log_likelihood(double * x,
                                     double * kappa,
                                     int num_classes):
    '''
    Return log likelihood of dirichlet distribution.
    '''
    cdef int i
    cdef double log_likelihood
    
    log_likelihood = 0
        
    for i in range(num_classes):
        log_likelihood += (kappa[i] - 1) * log(x[i])
    
    return log_likelihood

#=======================================================================================================================
# General mixture model code.
#=======================================================================================================================
cdef double * mixture_posterior(double * log_likelihood,
                                double * mix_weight,
                                int num_classes):
    '''
    Compute normalised posterior probabilities from mixture distribution by adding log_likelihood and 
    then normalising.
    
    Allocates posterior which will need to be freed by caller.
    ''' 
    cdef int i
    cdef double * posterior = < double *> malloc(num_classes * sizeof(double))

    for i in range(num_classes):
        log_likelihood[i] = log(mix_weight[i]) + log_likelihood[i]
    
    log_space_normalise_row(log_likelihood, num_classes)
    
    for i in range(num_classes):
        posterior[i] = exp(log_likelihood[i])
    
    return posterior

#=======================================================================================================================
# Code for doing log space normalisation
#=======================================================================================================================
cdef void log_space_normalise_row(double * log_X, int size):
    '''
    Normalise log_X so that 
    
    exp(log_X[0]) + ... + exp(log_X[1]) == 1
    
    Done in place so log_X is modified.
    '''
    cdef int i
    cdef double norm_const
    
    norm_const = log_sum_exp(log_X, size)
    
    for i in range(size):
        log_X[i] = log_X[i] - norm_const    

cdef double log_sum_exp(double * log_X, int size):
    '''
    Given a c-array log_X of values compute log( exp(log_X[0]) + ... + exp(log_X[size]) ).
    
    Numerically safer than naive method.
    '''
    cdef int i
    cdef double max_exp, total
 
    max_exp = log_X[0]
 
    for i in range(size):
        if max_exp < log_X[i]:
            max_exp = log_X[i]

    total = 0
    for i in range(size):
        total += exp(log_X[i] - max_exp)
    
    return log(total) + max_exp
