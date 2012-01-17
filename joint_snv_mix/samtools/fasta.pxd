'''
Created on 2012-01-17

@author: Andrew Roth
'''
from libc.stdlib cimport free
from libc.string cimport strdup

from joint_snv_mix.samtools.samtools_clib cimport fai_load, fai_destroy, faidx_fetch_seq, faidx_t

cdef extern from "ctype.h":
    int toupper(int c)

cdef class FastaFile:
    cdef char * _file_name
    cdef faidx_t * _fasta_file
   
    cdef faidx_t * get_file_pointer(self)
    cdef char * get_reference_base(self, char * reference, int position)
    cdef char * _fetch(self, char * reference, int start, int end, int * length)
