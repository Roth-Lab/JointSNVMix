#cython: cdivision=True

from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow

from joint_snv_mix.classifiers.classifier cimport Classifier

cdef class ThresholdClassifier(Classifier):    
    cdef float _normal_threshold
    cdef float _tumour_threshold
    
    cdef int _get_genotype(self, int a, int b, float freq_threshold)