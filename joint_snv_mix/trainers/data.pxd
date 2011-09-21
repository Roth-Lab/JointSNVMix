cdef class SnvMixData(object):
    pass
 
cdef class SnvMixOneData(SnvMixData):
    cdef int counts[2]
 
cdef class SnvMixTwoData(SnvMixData):
    cdef int depth
    cdef double * q
    cdef double * r