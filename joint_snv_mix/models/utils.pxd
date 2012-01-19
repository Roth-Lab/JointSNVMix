'''
Created on 2011-08-04

@author: Andrew Roth
'''
from libc.math cimport exp, log

cdef double binomial_log_likelihood(int a, int b, double mu)
cdef double beta_log_likelihood(double mu, double a, double b)
cdef double dirichlet_log_likelihood(list x, list kappa)

cdef double snv_mix_two_log_likelihood(double * q, double * r, int d, double mu)
cdef double snv_mix_two_single_read_log_likelihood(double q, double r, double mu)
cdef snv_mix_two_expected_a(double q, double r, double mu)
cdef snv_mix_two_expected_b(double q, double r, double mu)

cdef void log_space_normalise(double * log_X, int size)
cdef double log_sum_exp(double * log_X, int size)