from libc.math cimport exp, log
from libc.stdlib cimport malloc, free

cdef double multinomial_log_likelihood(int * x, double * mu, int num_classes)

cdef double dirichlet_log_likelihood(double * x, double * kappa, int num_classes)

cdef double * mixture_posterior(double * log_likelihood, double * mix_weight, int num_classes)

cdef void log_space_normalise_row(double * log_X, int size)

cdef double log_sum_exp(double * log_X, int size)