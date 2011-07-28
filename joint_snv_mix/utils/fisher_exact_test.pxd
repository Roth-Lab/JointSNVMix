# This code was originally from https://github.com/brentp/fishers_exact_test.

#cython: cdivision=True
from joint_snv_mix.utils.special_functions cimport hypergeometric_probability

cdef extern from "math.h":
    double log(double)
    double exp(double)

cdef inline int imin2(int a, int b):
    return a if a < b else b

cdef inline int imax2(int a, int b):
    return a if a > b else b

cdef class PValues(object):
    cdef double left_tail
    cdef double right_tail
    cdef double two_tail

cpdef PValues fisher_exact_test(int a_true, int a_false, int b_true, int b_false)
