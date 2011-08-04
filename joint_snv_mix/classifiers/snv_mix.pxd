#cython: cdivision=True
from libc.stdlib cimport free
from libc.math cimport log

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.classifiers.shared cimport multinomial_log_likelihood, get_mixture_posterior, \
                                              combine_independent_probs

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMixClassifier(Classifier):    
    pass
             
cdef class SnvMixClassifierRefIterator(ClassifierRefIterator):   
    cdef double _log_mu_N[NUM_GENOTYPES][2]
    cdef double _log_mu_T[NUM_GENOTYPES][2]
    cdef double _log_pi_N[NUM_GENOTYPES]
    cdef double _log_pi_T[NUM_GENOTYPES]