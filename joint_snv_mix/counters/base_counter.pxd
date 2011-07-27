from csamtools cimport Samfile, Fastafile, IteratorColumnRegion, PileupProxy,\
                       bam1_t, bam1_seq, bam1_qual, bam_pileup1_t, bam_dup1, bam_destroy1

from joint_snv_mix.counters.counter cimport Counter, CounterRefIterator, CounterRow
from joint_snv_mix.counters.shared cimport counts_struct, base_counts_struct, strcmp
    
cdef extern from "stdint.h":
    ctypedef int uint8_t

DEF ASCII_OFFSET = 33

cdef class BaseCounter(Counter):
    cdef Samfile _bam_file
    cdef int _min_base_qual
    cdef int _min_map_qual
    
cdef class BaseCounterRow(CounterRow):
    cdef counts_struct _counts
    
    cdef int get_counts(self, char * base)
    cdef base_counts_struct get_base_counts(self, char * base)

cdef class BaseCounterRefIterator(CounterRefIterator):
    cdef int _min_base_qual
    cdef int _min_map_qual
    cdef IteratorColumnRegion _pileup_iter
    
    cdef counts_struct _parse_pileup_column(self, PileupProxy pileup_column)