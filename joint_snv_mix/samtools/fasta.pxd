'''
Created on 2012-01-17

@author: Andrew Roth
'''
cdef class FastaFile:
    cdef char * _file_name
    cdef faidx_t * _fasta_file
   
    cdef faidx_t * get_file_pointer(self)
    cdef char * get_reference_base(self, char * reference, int position)
    cdef char * _fetch(self, char * reference, int start, int end, int * length)
