'''
Created on 2012-01-25

@author: Andrew Roth
'''
from joint_snv_mix.counter_row cimport JointBinaryCounterRow

cdef class RefIterator(object):
    cdef int _pos
    cdef char * _ref

    cdef JointBinaryCounterRow _current_row

    cdef cnext(self)
    cdef advance_position(self)
    cdef parse_current_position(self)
