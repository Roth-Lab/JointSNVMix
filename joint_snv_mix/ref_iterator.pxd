'''
Created on 2012-01-25

@author: Andrew Roth
'''
from joint_snv_mix.counter cimport JointBinaryCounterRow

cdef class RefIterator(object):
    cdef JointBinaryCounterRow _current_row

    cdef cnext(self)
    cdef advance_position(self)
    cdef parse_current_position(self)
