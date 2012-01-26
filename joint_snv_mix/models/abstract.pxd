'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

from libc.math cimport exp, log
from libc.stdlib cimport malloc, free

from joint_snv_mix.counter_row cimport JointBinaryData, JointBinaryCountData

from joint_snv_mix.models.utils cimport log_space_normalise, log_sum_exp 


cdef class Priors(object):
    cdef tuple _pi
 
cdef class Parameters(object):
    cdef tuple _pi

    cdef update(self, Ess ess, Priors priors)    
    cdef _update_pi(self, n, prior)   

cdef class MixtureModel(object):
    cdef int _num_joint_genotypes    
    cdef double * _resp
    
    cdef Density _density
    cdef Ess _ess
    
    cdef Priors _priors
    cdef Parameters _params

    cdef _predict(self, JointBinaryData data_point)
    
    cdef _E_step(self, data)
    cdef _M_step(self)

    cdef double _get_log_likelihood(self, data)
    cdef _get_prior_log_likelihood(self)
    
cdef class Density(object):  
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
     
    cdef get_responsibilities(self, JointBinaryData data_point, double * resp)        
    
    cdef double get_log_likelihood(self, JointBinaryData data_point, double * ll)
        
    cdef set_params(self, Parameters params)    

    cdef _get_complete_log_likelihood(self, JointBinaryData data_point, double * ll)

cdef class Ess(object):
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _n

    cdef reset(self)
    cdef set_params(self, Parameters params)       
    cdef update(self, JointBinaryData data_point, double * resp)
