'''
Created on 2012-01-22

@author: Andrew Roth
'''
from libc.stdio cimport fopen, fclose, fputs, fprintf, stdout, FILE
from libc.stdlib cimport free
from libc.string cimport strdup

from joint_snv_mix.counter_row cimport JointBinaryCounterRow, JointBinaryData

cdef class CResultsWriter(object):
    cdef FILE * _file_ptr
    cdef char * _format_string

    cdef _init_format_string(self)        
    
    cdef close(self)
    cdef write_position(self, JointBinaryCounterRow row, double * probs)
