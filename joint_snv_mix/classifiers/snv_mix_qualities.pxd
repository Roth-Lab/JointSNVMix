#cython: cdivision=True
from libc.stdlib cimport free
from libc.math cimport log

from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow, \
                                                                 JointBinaryQualityCounterIterator, \
                                                                 base_map_qualities_struct
                                                                 
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator

from joint_snv_mix.classifiers.shared cimport snv_mix_2_log_likelihood, get_mixture_posterior, \
                                              combine_independent_probs, get_phred_to_prob_qual_map

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMix2Classifier(Classifier):    
    pass
             
cdef class SnvMix2ClassifierRefIterator(ClassifierRefIterator):              
    cdef double * _qual_map
       
    cdef double _mu_N[NUM_GENOTYPES]
    cdef double _mu_T[NUM_GENOTYPES]
    cdef double _log_pi_N[NUM_GENOTYPES]
    cdef double _log_pi_T[NUM_GENOTYPES]