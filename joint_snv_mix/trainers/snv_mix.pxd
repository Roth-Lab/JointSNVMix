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

from joint_snv_mix.utils.log_pdf cimport dirichlet_log_likelihood, log_space_normalise_row, log_sum_exp

DEF NUM_GENOTYPES = 3
DEF NUM_BASES = 2
 
cdef class SnvMixCpt(object):
    cdef double * get_resp(self)
    cdef double * get_expected_counts_a(self)
    cdef double * get_expected_counts_b(self)
    cdef double get_log_sum(self)
 
cdef class SnvMixOneCpt(SnvMixCpt):
    cdef int _a
    cdef int _b
    cdef double * _cpt_array
    
    cdef _init_cpt_array(self, SnvMixOneData data, SnvMixParameters params)
    cdef double _binomial_log_likelihood(self, int a, int b, double mu)

      
cdef class SnvMixTwoCpt(SnvMixCpt):
    cdef int _depth
    cdef double **** _cpt_array
    cdef double _marginal
    cdef double _resp[NUM_GENOTYPES]

    cdef double * _get_expected_counts(self, int a)    
    cdef _init_cpt_array(self, SnvMixTwoData data, SnvMixParameters params)
    cdef _fill_cpt_array(self, SnvMixTwoData data, SnvMixParameters params)    
    cdef double _get_read_complete_likelihood(self, int a, int z, double q, double r, double mu)        
    cdef void _make_cpt_array(self)
    cdef void _free_cpt_array(self)
    cdef double ** _get_read_marginals(self, SnvMixTwoData data, SnvMixParameters params)
    cdef void _free_read_marginals(self, double ** read_marginals)
    cdef double * _get_class_marginals(self, double ** read_marginals, SnvMixParameters params)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixEss(object):
    cdef double _a[NUM_GENOTYPES]
    cdef double _b[NUM_GENOTYPES]
    cdef double _n[NUM_GENOTYPES]
   
    cdef void reset(self)
    cdef update(self, SnvMixCpt cpt)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixModel(object):
    cdef SnvMixParameters _params

    cdef double _get_lower_bound(self, list data)
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data)
    cdef double _get_log_likelihood(self, SnvMixData data)
    
cdef class PairedSnvMixModel(object):
    cdef PairedSnvMixParameters _params
    
    cdef SnvMixModel _normal_model
    cdef SnvMixModel _tumour_model
    
cdef class PairedSnvMixOneModel(PairedSnvMixModel):
    pass

cdef class PairedSnvMixTwoModel(PairedSnvMixModel):
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

