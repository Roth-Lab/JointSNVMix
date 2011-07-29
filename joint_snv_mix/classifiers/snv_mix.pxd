#cython: cdivision=True

from libc.math cimport exp, log

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test
from joint_snv_mix.utils.normalise cimport log_space_normalise_row

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMixClassifier(Classifier):    
    pass
             
cdef class SnvMixClassifierRefIterator(ClassifierRefIterator):   
    cdef double _log_mu_N[NUM_GENOTYPES][2]
    cdef double _log_mu_T[NUM_GENOTYPES][2]
    cdef double _log_pi_N[NUM_GENOTYPES]
    cdef double _log_pi_T[NUM_GENOTYPES]

    cdef void _get_probs(
                         self, 
                         double probs[NUM_GENOTYPES], 
                         int a, 
                         int b, 
                         double mu[NUM_GENOTYPES][2], 
                         double pi[NUM_GENOTYPES]
                         )
    
    cdef void _combine_probs(
                             self, 
                             double joint_probs[NUM_JOINT_GENOTYPES], 
                             double normal_probs[NUM_GENOTYPES], 
                             double tumour_probs[NUM_GENOTYPES]
                             )
