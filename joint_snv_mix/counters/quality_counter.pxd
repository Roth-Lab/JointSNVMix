from libc.stdlib cimport malloc, free

from csamtools cimport Samfile, Fastafile, IteratorColumnRegion, PileupProxy, \
                       bam1_t, bam1_seq, bam1_qual, bam_pileup1_t, bam_dup1, bam_destroy1

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport SingleSampleCounterRow
from joint_snv_mix.counters.ref_iterator cimport CRefIterator, RefIterator
from joint_snv_mix.counters.shared cimport base_counts_struct, column_struct
    
cdef extern from "stdint.h":
    ctypedef int uint8_t

DEF ASCII_OFFSET = 33
DEF NUM_QUAL_VAL = 256

cdef class QualityCounter(Counter):
    cdef Samfile _bam_file
    
cdef class QualityCounterRow(SingleSampleCounterRow):
    cdef char * _bases
    cdef int * _base_quals
    cdef int * _map_quals

cdef class QualityCounterRefIterator(RefIterator):
    cdef double _qual_map[NUM_QUAL_VAL]
    cdef CRefIterator _ref_iter

    cdef advance_position(self)
    cdef parse_current_position(self)
    
    
