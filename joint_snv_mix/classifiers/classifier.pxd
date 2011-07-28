from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryBaseCounter, JointBinaryBaseCounterIterator, JointBinaryCounterRow

cdef class Classifier(object):
    # List of references in Bam file(s) the counter works over.
    cdef tuple _refs
    
    cdef JointBinaryBaseCounter _counter

cdef class ClassifierRow(object):
    # Ref for the row
    cdef char * _ref
    
    # 0-based position of the row.
    cdef int _position
    
    cdef char * _ref_base
    
    cdef char * _non_ref_base
    
    cdef tuple _counts
    
    cdef tuple _labels  

cdef class ClassifierRefIterator(object):
    # Ref which the iterator runs over.
    cdef char * _ref
    
    # 0-based current position of the iterator.
    cdef int _position
    
    cdef object _current_row
    
    cdef int _min_normal_depth
    
    cdef int _min_tumour_depth
    
    cdef JointBinaryBaseCounterIterator _iter
    
    cdef cnext(self)
    
    cdef tuple _get_labels(self)    