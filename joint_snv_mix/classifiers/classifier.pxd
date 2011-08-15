from libc.stdlib cimport free
from libc.stdio cimport fopen, fclose, fputs, FILE, fprintf

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow

DEF NUM_JOINT_GENOTYPES = 9

cdef class ClassifierRow(object):
    # Ref for the row
    cdef char * _ref
    
    # 0-based position of the row.
    cdef int _position
    
    cdef char * _ref_base
    
    cdef char * _non_ref_base
    
    cdef int _counts[4]
    
    cdef double * _labels

cdef class Classifier(object):
    cdef ClassifierRow _classify(self, PairedSampleBinomialCounterRow row)    
    cdef double * _get_labels(self, PairedSampleBinomialCounterRow row)
