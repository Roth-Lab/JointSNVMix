'''
Created on 2012-01-25

@author: Andrew Roth
'''
from libc.stdio cimport fopen, fclose, fscanf, FILE, EOF, ftell, fseek
from libc.string cimport strcmp, strcpy

from joint_snv_mix.counter cimport JointBinaryCounterIterator
from joint_snv_mix.ref_iterator cimport RefIterator

cdef class PositionsCounter(object):    
    cdef char * _pos_file_name
    cdef dict index
    cdef tuple refs
       
    cdef _load_index(self)        
        
cdef class PositionsCounterRefIterator(RefIterator):
    cdef char * _ref
    cdef JointBinaryCounterIterator _ref_iter
    cdef PositionsIterator _pos_iter
        
    cdef advance_position(self)
    cdef parse_current_position(self)                 

cdef class PositionsIterator(object):
    cdef int _stop        
    cdef FILE * _file_p

    cdef cnext(self)