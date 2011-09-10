'''
Created on 2011-08-11

@author: Andrew Roth
'''
from libc.stdio cimport fopen, fclose, fscanf, FILE, EOF, ftell, fseek
from libc.string cimport strcmp, strcpy

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.ref_iterator cimport RefIterator

cdef class PositionsCounter(Counter):
    cdef dict _intervals
    cdef char * _pos_file_name
    cdef Counter _counter

    cdef _load_index(self)

cdef class PositionsIterator(object):
    cdef FILE * _file_p
    cdef int _stop
    cdef int _position
    
    cdef cnext(self)
    
cdef class PositionsCounterRefIterator(RefIterator):
    cdef RefIterator _ref_iter
    cdef PositionsIterator _pos_iter