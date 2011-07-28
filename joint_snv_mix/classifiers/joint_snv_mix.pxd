#cython: cdivision=True

from libc.math cimport exp, log

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test

cdef class JointSnvMixClassifier(Classifier):    
    pass
             
cdef class JointSnvMixClassifierRefIterator(ClassifierRefIterator):   
    cdef float _mu_N[3]
    cdef float _mu_T[3]
    cdef float _pi[9]

    cdef void _compute_likelihood(self, float likelihood[3], int a, int b, float mu[3])
    cdef void _compute_joint_probs(
                                   self,
                                   float joint_probs[9],
                                   float normal_likelihood[3],
                                   float tumour_likelihood[3],
                                   float pi[9]
                                   )
    cdef void _normalise_log_probs(self, float probs[9])
    cdef float _compute_log_norm_constant(self, float probs[9])
