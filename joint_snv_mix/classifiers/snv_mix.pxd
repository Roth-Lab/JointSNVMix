#cython: cdivision=True

from libc.math cimport exp, log

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test

cdef class SnvMixClassifier(Classifier):    
    pass
             
cdef class SnvMixClassifierRefIterator(ClassifierRefIterator):   
    cdef tuple _mu_N
    cdef tuple _mu_T
    cdef tuple _pi_N
    cdef tuple _pi_T

    cdef tuple _get_probs(self, int a, int b, tuple mu, tuple pi)
    cdef void _normalise_log_probs(self, float probs[3])
    cdef float _compute_log_norm_constant(self, float probs[3])
    cdef tuple _combine_probs(self, tuple normal_probs, tuple tumour_probs)
