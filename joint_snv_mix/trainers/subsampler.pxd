cdef class PairedDataSubSampler(object):
    cdef int _skip_size
    cdef int _min_normal_depth
    cdef int _min_tumour_depth

cdef class Subsample(object):
    cdef list _normal_rows
    cdef list _tumour_rows 

    cdef add_row(self, PairedSampleBinomialCounterRow row)
     
cdef class SnvMixOneSubsample(Subsample):
    pass
 
cdef class SnvMixTwoSubsample(Subsample):
    pass
