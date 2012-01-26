'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

from libc.math cimport exp, log
from libc.stdlib cimport malloc, free

from joint_snv_mix.counter_row cimport JointBinaryData, JointBinaryCountData

from joint_snv_mix.models.abstract cimport Density, Ess, MixtureModel, Parameters, Priors

from joint_snv_mix.models.utils cimport beta_binomial_log_likelihood, dirichlet_log_likelihood, \
                                        log_space_normalise, log_sum_exp 

cdef class BetaBinomialPriors(Priors):
    pass
 
cdef class BetaBinomialParameters(Parameters):
    cdef tuple _alpha_N
    cdef tuple _alpha_T
    
    cdef tuple _beta_N
    cdef tuple _beta_T

cdef class BetaBinomialModel(MixtureModel):   
    pass
    
cdef class BetaBinomialDensity(Density):
    cdef double * _alpha_N
    cdef double * _alpha_T
        
    cdef double * _beta_N
    cdef double * _beta_T
        
    cdef double * _log_mix_weights

    cdef _init_arrays(self)

cdef class BetaBinomialEss(Ess):
    pass
