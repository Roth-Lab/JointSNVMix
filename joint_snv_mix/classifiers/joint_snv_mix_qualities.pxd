DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

from libc.stdlib cimport free
from libc.math cimport log

from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow
                                                                 
from joint_snv_mix.classifiers.classifier cimport Classifier

from joint_snv_mix.classifiers.shared cimport snv_mix_2_log_likelihood, get_joint_posterior, \
                                              get_phred_to_prob_qual_map

cdef class JointSnvMixTwoClassifier(Classifier):
    cdef double * _qual_map
       
    cdef double _mu_N[NUM_GENOTYPES]
    cdef double _mu_T[NUM_GENOTYPES]
    cdef double _log_pi[NUM_JOINT_GENOTYPES]
