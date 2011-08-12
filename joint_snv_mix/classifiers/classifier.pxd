from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.ref_iterator cimport JointRefIterator 

cdef class Classifier(object):
    cdef ClassifierRow _classify(self, PairedSampleBinomialCounterRow row)
    cdef tuple _get_labels(self, PairedSampleBinomialCounterRow row)

cdef class ClassifierRefIterator(object):
    # Ref which the iterator runs over.
    cdef char * _ref
    
    # 0-based current position of the iterator.
    cdef int _position
    
    cdef object _current_row
    
    cdef int _min_normal_depth
    
    cdef int _min_tumour_depth
    
    cdef JointRefIterator _iter
    
    cdef cnext(self)
    
    cdef tuple _get_labels(self)    
