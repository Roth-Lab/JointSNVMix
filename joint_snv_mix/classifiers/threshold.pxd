#cython: cdivision=True

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

cdef class ThresholdClassifier(Classifier):    
    pass
             
cdef class ThresholdClassifierRefIterator(ClassifierRefIterator):   
    cdef float _normal_threshold
    cdef float _tumour_threshold
    
    cdef int _get_genotype(self, int a, int b, float freq_threshold)