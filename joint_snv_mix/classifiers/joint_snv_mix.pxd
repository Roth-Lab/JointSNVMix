#cython: cdivision=True
from libc.stdlib cimport free
from libc.math cimport log

from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow

from joint_snv_mix.classifiers.classifier cimport Classifier

from joint_snv_mix.classifiers.shared cimport multinomial_log_likelihood, get_joint_posterior

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class JointSnvMixOneClassifier(Classifier):  
    cdef double _log_mu_N[NUM_GENOTYPES][2]
    cdef double _log_mu_T[NUM_GENOTYPES][2]
    cdef double _log_pi[9]