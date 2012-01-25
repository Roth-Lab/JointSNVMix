'''
Created on 2011-08-04

@author: Andrew Roth
'''
#=======================================================================================================================
# Log likelihoods
#=======================================================================================================================
cpdef double binomial_log_likelihood(int a, int b, double mu):
    return a * log(mu) + b * log(1 - mu)

cpdef double beta_binomial_log_likelihood(int a, int b, double alpha, double beta):
    return log_beta(a + alpha, b + beta) - log_beta(alpha, beta)

cpdef double beta_log_likelihood(double mu, double a, double b):
    return (a - 1) * log(mu) + (b - 1) * log(1 - mu)

cpdef double dirichlet_log_likelihood(tuple x, tuple kappa):
    cdef double k_i, x_i, log_likelihood
    
    log_likelihood = 0
        
    for x_i, k_i in zip(x, kappa):
        log_likelihood += (k_i - 1) * log(x_i)
    
    return log_likelihood

cdef double snv_mix_two_log_likelihood(double * q, double * r, int d, double mu):
    cdef int i
    cdef double log_likelihood
    
    log_likelihood = 0
    
    for i in range(d):
        log_likelihood += log(snv_mix_two_single_read_likelihood(q[i], r[i], mu))
    
    return log_likelihood

cpdef double snv_mix_two_single_read_likelihood(double q, double r, double mu):
    return 0.5 * (1 - r) + r * ((1 - q) * (1 - mu) + q * mu)

cpdef double snv_mix_two_expected_a(double q, double r, double mu):
    cdef double numerator, denominator
    
    numerator = mu * (0.5 + r * (q - 0.5))
    denominator = snv_mix_two_single_read_likelihood(q, r, mu)
    
    return numerator / denominator

cpdef double snv_mix_two_expected_b(double q, double r, double mu):
    cdef double numerator, denominator
    
    numerator = (1 - mu) * (0.5 + r * (0.5 - q))
    denominator = snv_mix_two_single_read_likelihood(q, r, mu)
    
    return numerator / denominator

#=======================================================================================================================
# Code for doing log space normalisation
#=======================================================================================================================   
cdef void log_space_normalise(double * log_X, int size):
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

#=======================================================================================================================
# Special functions
#=======================================================================================================================
cdef double log_beta(double a, double b):
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    
