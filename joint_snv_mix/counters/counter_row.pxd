cdef class CounterRow(object):
    cdef char * _ref    
    cdef int _position

cdef class SingleSampleCounterRow(CounterRow):    
    cdef int _depth

cdef class PairedSampleCounterRow(object):
    cdef int _normal_depth
    cdef int _tumour_depth