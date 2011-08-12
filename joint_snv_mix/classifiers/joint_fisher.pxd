#cython: cdivision=True

from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow

from joint_snv_mix.classifiers.classifier cimport Classifier
from joint_snv_mix.classifiers.independent_fisher cimport IndependentFisherClassifier

from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test

cdef extern from "math.h":
    float floor(float x)

cdef class JointFisherClassifier(Classifier):    
    cdef float _p_value_threshold
    
    cdef IndependentFisherClassifier _indep_classifier
