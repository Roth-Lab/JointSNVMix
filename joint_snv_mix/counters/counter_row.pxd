from joint_snv_mix.counters.shared cimport binary_counts_struct, base_counts_struct, strcmp

cdef extern from * :
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdlib.h":
    void qsort (void * ARRAY, size_t COUNT, size_t SIZE, int (*COMPARE)(const_void * , const_void *))

cdef class CounterRow(object):
    cdef char * _ref    
    cdef int _position

cdef class SingleSampleCounterRow(CounterRow):    
    cdef int _depth
    
    cdef int get_counts(self, char * base)

cdef class PairedSampleBinomialCounterRow(CounterRow):
    cdef int _normal_depth
    cdef int _tumour_depth
    
    cdef char * _ref_base
    cdef char * _non_ref_base
    
cdef char * get_non_ref_base(char * ref_base, SingleSampleCounterRow normal_row, SingleSampleCounterRow tumour_row)
cdef binary_counts_struct get_binary_counts(char * ref_base, char * non_ref_base, SingleSampleCounterRow row)
