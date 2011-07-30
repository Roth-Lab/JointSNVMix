cdef class Counter(object):
    # List of references in Bam file(s) the counter works over.
    cdef tuple _refs