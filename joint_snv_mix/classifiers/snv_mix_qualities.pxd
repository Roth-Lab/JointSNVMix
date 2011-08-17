#cython: cdivision=True
from libc.stdlib cimport free
from libc.math cimport log

from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow
                                                                 
from joint_snv_mix.classifiers.classifier cimport Classifier
                                              
from joint_snv_mix.trainers.snv_mix cimport PairedSnvMixParameters, SnvMixTwoData, SnvMixTwoCpt, makeSnvMixTwoData

from joint_snv_mix.classifiers.shared cimport combine_independent_probs                                             

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMixTwoClassifier(Classifier):           
    cdef PairedSnvMixParameters _params
       
