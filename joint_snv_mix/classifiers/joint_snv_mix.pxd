#cython: cdivision=True
from libc.stdlib cimport free
from libc.math cimport log

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.classifiers.shared cimport multinomial_log_likelihood, get_joint_posterior

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class JointSnvMixClassifier(Classifier):    
    pass
             
cdef class JointSnvMixClassifierRefIterator(ClassifierRefIterator):   
    cdef double _log_mu_N[NUM_GENOTYPES][2]
    cdef double _log_mu_T[NUM_GENOTYPES][2]
    cdef double _log_pi[9]