'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

from libc.math cimport exp, log
from libc.stdlib cimport malloc, free

from joint_snv_mix.counter cimport JointBinaryData, JointBinaryCountData, JointBinaryQualityData
from joint_snv_mix.models.utils cimport binomial_log_likelihood, beta_log_likelihood, dirichlet_log_likelihood, \
                                        snv_mix_two_log_likelihood, snv_mix_two_expected_a, snv_mix_two_expected_b, \
                                        log_space_normalise, log_sum_exp 

cdef class JointSnvMixPriors(object):
    cdef tuple _mu_N
    cdef tuple _mu_T
    cdef tuple _pi
 
cdef class JointSnvMixParameters(object):
    cdef tuple _mu_N
    cdef tuple _mu_T
    cdef tuple _pi

cdef class JointSnvMixModel(object):
    cdef JointSnvMixPriors _priors
    cdef JointSnvMixParameters _params
    
    cdef _JointSnvMixDensity _density
    cdef _JointSnvMixEss _ess
    
    cdef int _num_joint_genotypes
    cdef double * _resp

    cdef _predict(self, JointBinaryData data_point)
    
    cdef _E_step(self, data)
    cdef _M_step(self)

    cdef _get_updated_mu(self, a, b, prior)            
    cdef _get_updated_pi(self, n, prior)

    cdef double _get_log_likelihood(self, data)
    cdef _get_prior_log_likelihood(self)

cdef class _JointSnvMixDensity(object):
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _mu_N
    cdef double * _mu_T
    cdef double * _log_mix_weights

    cdef _get_complete_log_likelihood(self, JointBinaryData data_point, double * ll)

    cdef _init_arrays(self)    
    
    cdef get_responsibilities(self, JointBinaryData data_point, double * resp)
    cdef double get_log_likelihood(self, JointBinaryData data_point, double * ll)    
    cdef set_params(self, JointSnvMixParameters params)
     
cdef class _JointSnvMixOneDensity(_JointSnvMixDensity):
    pass
                
 
cdef class _JointSnvMixTwoDensity(_JointSnvMixDensity):
    pass

cdef class _JointSnvMixEss(object):
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _a_N
    cdef double * _b_N
    cdef double * _a_T
    cdef double * _b_T
    cdef double * _n
    
    cdef double * _mu_N
    cdef double * _mu_T

    cdef update(self, JointBinaryData data_point, double * resp)

    cdef _init_arrays(self)
    cdef reset(self)
    cdef set_params(self, JointSnvMixParameters params)
 
cdef class _JointSnvMixOneEss(_JointSnvMixEss):
    pass    
 
cdef class _JointSnvMixTwoEss(_JointSnvMixEss):
    pass