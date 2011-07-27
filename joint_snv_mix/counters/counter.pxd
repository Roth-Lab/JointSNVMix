cdef class Counter(object):
    # List of references in Bam file(s) the counter works over.
    cdef tuple _refs

cdef class CounterRow(object):
    # Ref for the row
    cdef char * _ref
    
    # 0-based position of the row.
    cdef int _position    

cdef class CounterRefIterator(object):
    # Ref which the iterator runs over.
    cdef char * _ref
    
    # 0-based current position of the iterator.
    cdef int _position
    
    cdef object _current_row
    
    cdef cnext(self)