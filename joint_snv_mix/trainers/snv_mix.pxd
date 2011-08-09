'''
Created on 2011-08-04

@author: Andrew Roth
'''
from libc.stdlib cimport free, malloc
from libc.math cimport abs, exp, log

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow
from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow
from joint_snv_mix.counters.ref_iterator cimport RefIterator

from joint_snv_mix.counters.shared cimport binary_counts_struct, base_map_qualities_struct

from joint_snv_mix.trainers.trainer cimport Trainer
from joint_snv_mix.utils.log_pdf cimport multinomial_log_likelihood, dirichlet_log_likelihood, log_space_normalise_row

DEF NUM_GENOTYPES = 3
DEF NUM_BASES = 2

cdef class SnvMixData(object):
    pass
 
cdef class SnvMixOneData(SnvMixData):
    cdef int counts[2]
 
cdef class SnvMixTwoData(SnvMixData):
    cdef int depth
    cdef int * labels
    cdef double * q
    cdef double * r    

#---------------------------------------------------------------------------------------------------------------------- 
cdef class PairedDataSubSampler(object):
    cdef int _skip_size
    cdef int _min_normal_depth
    cdef int _min_tumour_depth
    
    cdef _add_row_to_sample(self, dict sample, PairedSampleBinomialCounterRow row)
     
cdef class SnvMixOneSubsampler(PairedDataSubSampler):
    pass
 
cdef class SnvMixTwoSubsampler(PairedDataSubSampler):
    pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixPriors(object):
    cdef double mu[NUM_GENOTYPES][2]
    cdef double pi[NUM_GENOTYPES]
 
cdef class PairedSnvMixPriors(object):
    cdef SnvMixPriors _normal_priors
    cdef SnvMixPriors _tumour_priors

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixParameters(object):
    cdef SnvMixPriors priors
    
    cdef double mu[NUM_GENOTYPES]
    cdef double pi[NUM_GENOTYPES]

    cdef update(self, double * n, double * a, double *)
    
    cdef _update_mu(self, double * a, double * b)            
    cdef _update_pi(self, double * n)
    cdef double _get_prior_log_likelihood(self)
 
cdef class PairedSnvMixParameters(object):
    cdef SnvMixParameters _normal_params
    cdef SnvMixParameters _tumour_params

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixCpt(object):
    cdef double * get_resp(self)
    cdef double * get_expected_counts_a(self)
    cdef double * get_expected_counts_b(self)
    cdef double marginalise(self)
 
cdef class SnvMixOneCpt(SnvMixCpt):
    cdef int _a
    cdef int _b
    cdef double * _cpt_array
    
    cdef _init_cpt_array(self, SnvMixOneData data, SnvMixParameters params)
    cdef double _binomial_log_likelihood(self, int a, int b, double mu)

      
cdef class SnvMixTwoCpt(SnvMixCpt):
    cdef int _depth
    cdef double **** _cpt_array

    cdef double * _get_expected_counts(self, int a)    
    cdef _init_cpt_array(self, SnvMixTwoData data, SnvMixParameters params)
    cdef _fill_cpt_array(self, SnvMixTwoData data, SnvMixParameters params)    
    cdef _normalise_cpt_array(self)
    cdef double _get_read_complete_likelihood(self, int a, int z, double q, double r, double mu)        
    cdef void _make_cpt_array(self)
    cdef void _free_cpt_array(self)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixEss(object):
    cdef double a[NUM_GENOTYPES]
    cdef double b[NUM_GENOTYPES]
    cdef double n[NUM_GENOTYPES]
   
    cdef void reset(self)
    cdef update(self, SnvMixCpt cpt)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixModel(object):
    cdef SnvMixParameters params

    cdef double _get_lower_bound(self, list data)
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data)
    cdef double _get_log_likelihood(self, SnvMixData data)
    
cdef class PairedSnvMixModel(object):
    pass
 
cdef class SnvMixOneModel(SnvMixModel):
    pass

cdef class SnvMixTwoModel(SnvMixModel):
    pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixModelTrainer(object):
    cdef bint _converged
    cdef int _iters
    cdef double _convergence_threshold
    cdef int _max_iters

    cdef list _lower_bounds
    cdef SnvMixModel _model

    cdef train(self, list data)
    cdef SnvMixEss _do_e_step(self, list data)
    cdef void _do_m_step(self, SnvMixEss ess)
    cdef _check_convergence(self, list data)

