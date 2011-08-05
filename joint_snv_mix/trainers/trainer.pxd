'''
Created on 2011-08-04

@author: Andrew Roth
'''
from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.ref_iterator cimport RefIterator

cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng_type:
        pass
    ctypedef struct gsl_rng:
        pass
    
    gsl_rng_type * gsl_rng_mt19937
    
    gsl_rng * gsl_rng_alloc(gsl_rng_type * T)
    
    double gsl_rng_uniform (gsl_rng * r)

cdef class Trainer(object):
    cdef list _load_data_set(self, Counter counter, double sub_sample_fraction)
