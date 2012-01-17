'''
Created on 2012-01-17

@author: Andrew Roth
'''
cdef class BamFile:
    cdef char * _file_name
    cdef samfile_t * _bam_file
    cdef bam_index_t * _index
    
    cdef samfile_t * get_file_pointer(self)
    cdef bam_index_t * get_index(self)