'''
Classes for creating iterators for count data over a pair of genomes. 

Created on 2012-01-18

@author: Andrew Roth
'''
from libc.stdlib cimport malloc, free
from libc.string cimport strcmp, strdup

from joint_snv_mix.samtools.bam cimport BamFile
from joint_snv_mix.samtools.fasta cimport FastaFile
from joint_snv_mix.samtools.pileup cimport PileupColumn, PileupIterator

cdef extern from * :
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdlib.h":
    void qsort (void * ARRAY, size_t COUNT, size_t SIZE, int (*COMPARE)(const_void * , const_void *))

ctypedef struct base_counts_struct:
    char * base
    int counts

cdef class JointBinaryCounter(object):
    cdef char * _type
    
    cdef BamFile _normal_bam
    cdef BamFile _tumour_bam
    
    cdef FastaFile _ref_genome
    
    cdef tuple _refs 
        
cdef class JointBinaryCounterIterator(object):
    cdef char * _type
    cdef char * _ref

    cdef int _min_base_qual
    cdef int _min_map_qual
    
    cdef PileupIterator _normal_iter
    cdef PileupIterator _tumour_iter
    
    cdef FastaFile _ref_genome
    
    cdef int _pos
    
    cdef JointBinaryCounterRow _current_row
    
    cdef cnext(self)
    cdef advance_position(self)
    cdef parse_current_position(self)
    
    cdef JointBinaryCounterRow _make_counter_row(self, PileupColumn normal_column, PileupColumn tumour_column)
    
    cdef JointBinaryData _make_count_data(self, char * ref_base, char * var_base,
                          PileupColumn normal_column, PileupColumn tumour_column)
    
    cdef JointBinaryData _make_quality_data(self, char * ref_base, char * var_base,
                                            PileupColumn normal_column, PileupColumn tumour_column)

cdef class JointBinaryCounterRow(object):
    cdef char * _ref
    cdef int _pos
    
    cdef char * _ref_base
    cdef char * _var_base
    
    cdef JointBinaryData _data
    
cdef class JointBinaryBaseCounterRow(JointBinaryCounterRow):
    pass

cdef class JointBinaryQualityCounterRow(JointBinaryCounterRow):
    pass

#=======================================================================================================================
# Data object 
#=======================================================================================================================
cdef class JointBinaryData(object):
    cdef int _a_N
    cdef int _b_N
    
    cdef int _a_T
    cdef int _b_T

cdef class JointBinaryCountData(JointBinaryData):
    pass

cdef class JointBinaryQualityData(JointBinaryData):
    cdef int _d_N
    cdef int _d_T

    cdef double * _q_N
    cdef double * _r_N
    
    cdef double * _q_T
    cdef double * _r_T
