'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

from libc.math cimport exp, log
from libc.stdlib cimport malloc, free

from joint_snv_mix.counter cimport JointBinaryData, JointBinaryQualityData

from joint_snv_mix.models.abstract cimport Density, Ess, MixtureModel
from joint_snv_mix.models.binomial cimport BinomialParameters, BinomialPriors

from joint_snv_mix.models.utils cimport beta_log_likelihood, dirichlet_log_likelihood, \
                                        snv_mix_two_log_likelihood, snv_mix_two_expected_a, snv_mix_two_expected_b, \
                                        log_space_normalise, log_sum_exp 

cdef class SnvMixTwoModel(MixtureModel):
    cdef _get_updated_mu(self, a, b, prior)            

cdef class SnvMixTwoDensity(Density):
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _mu_N
    cdef double * _mu_T
    cdef double * _log_mix_weights

    cdef _init_arrays(self)

cdef class SnvMixTwoEss(Ess):
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

    cdef _init_arrays(self)
