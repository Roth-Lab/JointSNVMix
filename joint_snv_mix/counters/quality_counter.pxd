from libc.stdlib cimport malloc, free

from csamtools cimport Samfile, Fastafile, IteratorColumnRegion, PileupProxy,\
                       bam1_t, bam1_seq, bam1_qual, bam_pileup1_t, bam_dup1, bam_destroy1

from joint_snv_mix.counters.counter cimport Counter, CounterRefIterator, CounterRow
from joint_snv_mix.counters.shared cimport base_counts_struct
    
cdef extern from "stdint.h":
    ctypedef int uint8_t

DEF ASCII_OFFSET = 33
DEF NUM_QUAL_VAL = 256

cdef class QualityCounter(Counter):
    cdef Samfile _bam_file
    
cdef class QualityCounterRow(CounterRow):
    cdef int _num_reads
    cdef char * _bases
    cdef double * _base_quals
    cdef double * _map_quals
    
    cdef base_counts_struct get_base_counts(self, char * base)
    cdef int get_counts(self, char * base)

cdef class QualityCounterRefIterator(CounterRefIterator):
    cdef double _qual_map[NUM_QUAL_VAL]
    cdef IteratorColumnRegion _pileup_iter
    
    cdef _init_qual_map(self)
    cdef QualityCounterRow _parse_pileup_column(self, PileupProxy pileup_column)    
    
    