from libc.stdlib cimport free

from csamtools cimport Samfile, Fastafile

from joint_snv_mix.counters.counter cimport Counter, CounterRefIterator, CounterRow
from joint_snv_mix.counters.base_counter cimport BaseCounter, BaseCounterRow, BaseCounterRefIterator
from joint_snv_mix.counters.shared cimport binary_counts_struct, base_counts_struct, strcmp

cdef extern from * :
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdlib.h":
    void qsort (void * ARRAY, size_t COUNT, size_t SIZE, int (*COMPARE)(const_void * , const_void *))
    
cdef class JointBinaryBaseCounter(Counter):
    cdef BaseCounter _normal_counter
    cdef BaseCounter _tumour_counter    
    cdef Fastafile _ref_genome_fasta    
    
cdef class JointBinaryCounterRow(CounterRow):
    cdef char * _ref_base
    cdef char * _non_ref_base
    cdef binary_counts_struct _normal_counts
    cdef binary_counts_struct _tumour_counts

cdef class JointBinaryBaseCounterIterator(CounterRefIterator):
    cdef BaseCounterRefIterator _normal_iter
    cdef BaseCounterRefIterator _tumour_iter   
    cdef Fastafile _ref_genome_fasta    
    
    cdef _set_current_row(self)