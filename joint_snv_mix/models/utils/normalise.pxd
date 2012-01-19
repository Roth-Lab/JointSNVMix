'''
Created on 2011-07-28

@author: Andrew Roth
'''
from libc.math cimport exp, log

cdef void log_space_normalise_row(double * log_X, int size)

cdef double log_sum_exp(double * log_X, int size)
