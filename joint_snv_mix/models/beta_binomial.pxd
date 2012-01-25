'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

from libc.math cimport exp, log
from libc.stdlib cimport malloc, free

from joint_snv_mix.counter cimport JointBinaryData, JointBinaryCountData
from joint_snv_mix.models.utils cimport beta_binomial_log_likelihood, dirichlet_log_likelihood, \
                                        log_space_normalise, log_sum_exp 

cdef class BetaBinomialPriors(object):
    cdef tuple _pi
 
cdef class BetaBinomialParameters(object):
    cdef tuple _alpha_N
    cdef tuple _alpha_T
    
    cdef tuple _beta_N
    cdef tuple _beta_T
    
    cdef tuple _pi

cdef class BetaBinomialModel(object):
    cdef BetaBinomialPriors _priors
    cdef BetaBinomialParameters _params
    
    cdef _Density _density
    cdef _Ess _ess
    
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _resp

    cdef _predict(self, JointBinaryData data_point)
    
    cdef _E_step(self, data)
    cdef _M_step(self)

    cdef _get_updated_pi(self, n, prior)

    cdef double _get_log_likelihood(self, data)
    cdef _get_prior_log_likelihood(self)
    
cdef class _Density(object):
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _alpha_N
    cdef double * _alpha_T
    cdef double * _beta_N
    cdef double * _beta_T
    cdef double * _log_mix_weights
        
    cdef get_responsibilities(self, JointBinaryData data_point, double * resp)    
    
    cdef double get_log_likelihood(self, JointBinaryData data_point, double * ll)            
    
    cdef set_params(self, BetaBinomialParameters params)    

    cdef _init_arrays(self)
    cdef _get_complete_log_likelihood(self, JointBinaryData data_point, double * ll)

cdef class _Ess(object):
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes    
    cdef double * _n
    
    cdef update(self, JointBinaryData data_point, double * resp)    
    cdef reset(self)
    
