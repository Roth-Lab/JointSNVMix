from libc.stdlib cimport malloc, free

from csamtools cimport Samfile, Fastafile

from joint_snv_mix.counters.counter cimport Counter, CounterRefIterator, CounterRow
from joint_snv_mix.counters.quality_counter cimport QualityCounter, QualityCounterRow, QualityCounterRefIterator
from joint_snv_mix.counters.shared cimport binary_counts_struct, base_counts_struct, strcmp

cdef extern from * :
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdlib.h":
    void qsort (void * ARRAY, size_t COUNT, size_t SIZE, int (*COMPARE)(const_void * , const_void *))

ctypedef struct binary_depth_struct:
    int A
    int B

ctypedef struct binary_quality_struct:
    double * A
    double * B

ctypedef struct base_map_qualities_struct:
    binary_depth_struct depth
    binary_quality_struct base_quals
    binary_quality_struct map_quals
    
cdef class JointBinaryQualityCounter(Counter):
    cdef QualityCounter _normal_counter
    cdef QualityCounter _tumour_counter    
    cdef Fastafile _ref_genome_fasta    
    
cdef class JointBinaryQualityCounterRow(CounterRow):
    cdef char * _ref_base
    cdef char * _non_ref_base
    cdef base_map_qualities_struct _normal_data
    cdef base_map_qualities_struct _tumour_data

cdef class JointBinaryQualityCounterIterator(CounterRefIterator):
    cdef QualityCounterRefIterator _normal_iter
    cdef QualityCounterRefIterator _tumour_iter   
    cdef Fastafile _ref_genome_fasta    
    cdef QualityCounterRow _normal_row
    cdef QualityCounterRow _tumour_row
    
    cdef _set_current_row(self)
