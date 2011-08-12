#cython: cdivision=True

from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow, JointBinaryBaseCounterIterator
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow
from joint_snv_mix.classifiers.independent_fisher cimport IndependentFisherClassifierRefIterator
from joint_snv_mix.utils.fisher_exact_test cimport PValues, fisher_exact_test

cdef extern from "math.h":
    float floor(float x)

cdef class JointFisherClassifier(Classifier):    
    pass
             
cdef class JointFisherClassifierRefIterator(ClassifierRefIterator):
    cdef float _p_value_threshold
    
    cdef IndependentFisherClassifierRefIterator _indep_iter
