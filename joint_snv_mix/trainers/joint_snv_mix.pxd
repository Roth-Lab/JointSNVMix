'''
Created on 2011-08-04

Code for jsm

@author: Andrew Roth
'''
from libc.stdlib cimport free, malloc
from libc.math cimport abs, exp, log

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.ref_iterator cimport RefIterator
from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow
from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow

from joint_snv_mix.utils.log_pdf cimport dirichlet_log_likelihood, log_space_normalise_row, log_sum_exp
from joint_snv_mix.utils.special_functions cimport lncombination

from joint_snv_mix.trainers.snv_mix cimport SnvMixOneData, SnvMixTwoData, makeSnvMixOneData, makeSnvMixTwoData

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9
DEF NUM_BASES = 2

cdef class JointSnvMixData(object):
    pass
 
cdef class JointSnvMixOneData(JointSnvMixData):
    cdef SnvMixOneData _normal
    cdef SnvMixOneData _tumour
 
cdef class JointSnvMixTwoData(JointSnvMixData):
    cdef SnvMixTwoData _normal
    cdef SnvMixTwoData _tumour

#----------------------------------------------------------------------------------------------------------------------
cdef class PairedDataSubSampler(object):
    cdef int _skip_size
    cdef int _min_normal_depth
    cdef int _min_tumour_depth
    
    cdef _add_row_to_sample(self, list sample, PairedSampleBinomialCounterRow row)
 
cdef class JointSnvMixOneSubsampler(PairedDataSubSampler):
    pass
 
cdef class JointSnvMixTwoSubsampler(PairedDataSubSampler):
    pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixPriors(object):
    cdef double _mu_N[NUM_GENOTYPES][2]
    cdef double _mu_T[NUM_GENOTYPES][2]
    cdef double _pi[NUM_JOINT_GENOTYPES]

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixParameters(object):
    cdef JointSnvMixPriors _priors
    
    cdef double _mu_N[NUM_GENOTYPES]
    cdef double _mu_T[NUM_GENOTYPES]
    cdef double _pi[NUM_JOINT_GENOTYPES]

    cdef update(self, double * n, double * a_N, double * a_T, double * b_N, double * b_T)
    
    cdef _normalise_pi(self)
    cdef _update_mu(self, double * mu, double mu_prior[NUM_GENOTYPES][2], double * a, double * b)            
    cdef _update_pi(self, double * n)
    cdef double _get_prior_log_likelihood(self)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixCpt(object):
    cdef double get_log_sum(self)
    cdef double * get_resp(self)
    cdef double * get_expected_counts_a_N(self)
    cdef double * get_expected_counts_a_T(self)
    cdef double * get_expected_counts_b_N(self)
    cdef double * get_expected_counts_b_T(self)
 
cdef class JointSnvMixOneCpt(JointSnvMixCpt):
    cdef int _a_N
    cdef int _a_T
    cdef int _b_N
    cdef int _b_T
    
    cdef double _cpt_array[NUM_JOINT_GENOTYPES]
    
    cdef _init_cpt_array(self, JointSnvMixOneData data, JointSnvMixParameters params)
    cdef double _binomial_log_likelihood(self, int a, int b, double mu)
    cdef double * _get_normal_marginal_resp(self)
    cdef double * _get_tumour_marginal_resp(self)
    cdef double * _get_expected_counts(self, int counts, double * marginal_resp)

cdef class SampleCpt(object):
    cdef int _depth
    cdef double **** _cpt_array
    cdef double ** _read_marginals
    cdef double * _class_marginals

    cdef double * get_expected_counts_a(self, double * norm_const)
    cdef double * get_expected_counts_b(self, double * norm_const)
    
    cdef double * _get_expected_counts(self, int a, double * norm_const)
    
    cdef void _init_cpt_array(self, SnvMixTwoData data, double * mu)
    cdef double _get_read_complete_likelihood(self, int a, int z, double q, double r, double mu)
    
    cdef void _init_read_marginals(self)
    
    cdef void _init_class_marginals(self)
    
    cdef void _allocate_cpt_array(self)
    cdef void _free_cpt_array(self)
    
    cdef void _allocate_read_marginals(self)
    cdef void _free_read_marginals(self)
      
cdef class JointSnvMixTwoCpt(JointSnvMixCpt):
    cdef double _marginal
    cdef double _resp[NUM_JOINT_GENOTYPES]
    cdef double * _normal_counts_a
    cdef double * _normal_counts_b
    cdef double * _tumour_counts_a
    cdef double * _tumour_counts_b
    
    cdef double * _get_joint_class_marginals(self, double * normal_marginals, double * tumour_marginals, double * pi)
    
    cdef void _init_marginal(self, double * joint_marginals)
    cdef void _init_resp(self, double * joint_marginals)
    cdef void _init_normal_expected_counts(self, SampleCpt cpt, double * joint_marginals)   
    cdef void _init_tumour_expected_counts(self, SampleCpt cpt, double * joint_marginals)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixEss(object):
    cdef double _a_N[NUM_GENOTYPES]
    cdef double _a_T[NUM_GENOTYPES]
    cdef double _b_N[NUM_GENOTYPES]
    cdef double _b_T[NUM_GENOTYPES]
    cdef double _n[NUM_JOINT_GENOTYPES]
   
    cdef void reset(self)
    
    cdef update(self, JointSnvMixCpt cpt)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixModel(object):
    cdef JointSnvMixParameters _params

    cdef double _get_lower_bound(self, list data)
    cdef JointSnvMixCpt _get_complete_log_likelihood(self, JointSnvMixData data)
    cdef double _get_log_likelihood(self, JointSnvMixData data)
 
cdef class JointSnvMixOneModel(JointSnvMixModel):
    pass

cdef class JointSnvMixTwoModel(JointSnvMixModel):
    pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixModelTrainer(object):
    cdef bint _converged
    cdef int _iters
    cdef double _convergence_threshold
    cdef int _max_iters

    cdef list _lower_bounds
    cdef JointSnvMixModel _model

    cdef train(self, list data)
    cdef JointSnvMixEss _do_e_step(self, list data)
    cdef void _do_m_step(self, JointSnvMixEss ess)
    cdef _check_convergence(self, list data)

cdef JointSnvMixOneData makeJointSnvMixOneData(JointBinaryCounterRow row)
cdef JointSnvMixTwoData makeJointSnvMixTwoData(JointBinaryQualityCounterRow row)
