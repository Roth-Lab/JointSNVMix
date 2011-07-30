cdef class Counter(object):
    # List of references in Bam file(s) the counter works over.
    cdef tuple _refs

cdef class CounterRow(object):
    # Ref for the row
    cdef char * _ref
    
    # 0-based position of the row.
    cdef int _position