#cython: cdivision=True

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow, ClassifierOptions
from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test

cdef extern from "math.h":
    float floor(float x)

cdef class IndependentFisherClassifierOptions(ClassifierOptions):
    cdef int min_depth    
    cdef float min_var_freq
    cdef float hom_var_freq
    cdef float p_value_threshold
    cdef float expected_error_rate    

cdef class IndependentFisherClassifier(Classifier):    
    pass
             
cdef class IndependentFisherClassifierRefIterator(ClassifierRefIterator):
    cdef float _get_significance(self, int a, int b)
    cdef int _get_genotype(self, int a, int b)
