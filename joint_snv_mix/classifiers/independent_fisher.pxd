#cython: cdivision=True

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test

cdef extern from "math.h":
    float floor(float x)

cdef class IndependentFisherClassifier(Classifier):    
    pass
             
cdef class IndependentFisherClassifierRefIterator(ClassifierRefIterator):   
    cdef float _min_var_freq
    cdef float _hom_var_freq
    cdef float _p_value_threshold
    cdef float _expected_error_rate 
    cdef int _min_var_depth

    cdef float _get_significance(self, int a, int b)
    cdef int _get_genotype(self, int a, int b)
