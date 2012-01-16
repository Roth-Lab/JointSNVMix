DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

from libc.stdlib cimport free
from libc.math cimport log

from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow
from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow
                                                                 
from joint_snv_mix.classifiers.classifier cimport Classifier

from joint_snv_mix.trainers.joint_snv_mix cimport JointSnvMixParameters, JointSnvMixTwoCpt, JointSnvMixTwoData, \
                                                  makeJointSnvMixTwoData

cdef class JointSnvMixTwoClassifier(Classifier):
    cdef JointSnvMixParameters _params
