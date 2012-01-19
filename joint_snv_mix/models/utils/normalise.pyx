'''
Created on 2011-07-28

@author: Andrew Roth
'''
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
