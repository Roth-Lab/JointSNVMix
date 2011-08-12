from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow

cdef class ClassifierRow(object):
    # Ref for the row
    cdef char * _ref
    
    # 0-based position of the row.
    cdef int _position
    
    cdef char * _ref_base
    
    cdef char * _non_ref_base
    
    cdef tuple _counts
    
    cdef tuple _labels 

cdef class Classifier(object):
    cdef ClassifierRow _classify(self, PairedSampleBinomialCounterRow row)    
    cdef tuple _get_labels(self, PairedSampleBinomialCounterRow row)