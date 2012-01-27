'''
Created on 2012-01-18

@author: Andrew Roth
'''
from libc.math cimport pow
from libc.stdlib cimport malloc, free
from libc.string cimport strcmp, strdup

from joint_snv_mix.samtools.fasta cimport FastaFile
from joint_snv_mix.samtools.pileup cimport PileupColumn

#=======================================================================================================================
# External imports
#=======================================================================================================================
cdef extern from * :
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdlib.h":
    void qsort (void * ARRAY, size_t COUNT, size_t SIZE, int (*COMPARE)(const_void * , const_void *))

ctypedef struct base_counts_struct:
    char * base
    int counts

#=======================================================================================================================
# Classes
#=======================================================================================================================
cdef class RowFactory(object):
    cdef int _min_base_qual
    cdef int _min_map_qual
    
    cdef FastaFile _ref_genome
    
    cdef DataFactory _data_factory
    
    cdef JointBinaryCounterRow get_row(self, char * ref, int pos,
                                       PileupColumn normal_column, PileupColumn tumour_column)

    cdef char * _get_var_base(self, char * ref_base, PileupColumn normal_column, PileupColumn tumour_column)    

cdef class DataFactory(object):
    cdef int _min_base_qual
    cdef int _min_map_qual
    
    cdef JointBinaryData get_data(self, char * ref_base, char * var_base,
                                  PileupColumn normal_column, PileupColumn tumour_column)
    
cdef class CountDataFactory(DataFactory):
    pass

cdef class QualityDataFactory(DataFactory):
    cdef _get_aligment_probabilities(self,
                                     char * ref_base,
                                     char * var_base,
                                     double * q,
                                     double * r,
                                     PileupColumn column)

cdef class JointBinaryCounterRow(object):
    cdef char * _ref
    cdef int _pos
    
    cdef char * _ref_base
    cdef char * _var_base
    
    cdef JointBinaryData _data

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