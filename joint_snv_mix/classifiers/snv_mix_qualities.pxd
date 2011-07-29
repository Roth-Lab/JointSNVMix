#cython: cdivision=True

from libc.math cimport exp, log

from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow, \
                                                                 JointBinaryQualityCounterIterator, \
                                                                 JointBinaryQualityCounter, \
                                                                 base_map_qualities_struct
                                                                 
from joint_snv_mix.classifiers.classifier cimport Classifier, ClassifierRefIterator, ClassifierRow

from joint_snv_mix.utils.normalise cimport log_space_normalise_row

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMixClassifier(object):    
    cdef JointBinaryQualityCounter _counter
    cdef tuple _refs
             
cdef class SnvMixClassifierRefIterator(object):
    # Ref which the iterator runs over.
    cdef char * _ref
    
    # 0-based current position of the iterator.
    cdef int _position
    
    cdef object _current_row
    
    cdef JointBinaryQualityCounterIterator _iter
       
    cdef double _mu_N[NUM_GENOTYPES]
    cdef double _mu_T[NUM_GENOTYPES]
    cdef double _log_pi_N[NUM_GENOTYPES]
    cdef double _log_pi_T[NUM_GENOTYPES]

    cdef void _get_probs(
                         self,
                         double probs[NUM_GENOTYPES],
                         base_map_qualities_struct data,
                         double mu[NUM_GENOTYPES],
                         double log_pi[NUM_GENOTYPES]
                         )
    
    cdef void _combine_probs(
                             self,
                             double joint_probs[NUM_JOINT_GENOTYPES],
                             double normal_probs[NUM_GENOTYPES],
                             double tumour_probs[NUM_GENOTYPES]
                             )
    
    cdef double _compute_single_base_log_prob(self, double q, double r, double m)
    cdef tuple _get_labels(self)
    cdef cnext(self)
