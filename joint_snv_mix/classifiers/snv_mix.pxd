#cython: cdivision=True

from libc.math cimport exp, log

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test

cdef class SnvMixClassifier(Classifier):    
    pass
             
cdef class SnvMixClassifierRefIterator(ClassifierRefIterator):   
    cdef float _mu_N[3]
    cdef float _mu_T[3]
    cdef float _pi_N[3]
    cdef float _pi_T[3]

    cdef void _get_probs(self, float[3]probs, int a, int b, float mu[3], float pi[3])
    cdef void _normalise_log_probs(self, float probs[3])
    cdef float _compute_log_norm_constant(self, float probs[3])
    cdef void _combine_probs(self, float[9] joint_probs, float[3] normal_probs, float[3] tumour_probs)
