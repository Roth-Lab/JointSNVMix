'''
Created on 2012-01-26

@author: Andrew Roth
'''
from libc.math cimport pow, sqrt

from joint_snv_mix.samtools.bam cimport BamFile
from joint_snv_mix.models.binomial cimport BinomialPriors, BinomialParameters
from joint_snv_mix.models.snv_mix_two cimport SnvMixTwoModel
from joint_snv_mix.samtools.pileup cimport ExtendedPileupColumn, PileupIterator

from joint_snv_mix.counter_row cimport JointBinaryCounterRow

cdef extern from "math.h":
    double INFINITY

cdef class FeatureExtractor(object):
    cdef BamFile _normal_bam
    cdef BamFile _tumour_bam
    
    cdef SnvMixTwoModel _model
   
    cdef tuple _get_genome_features(self, JointBinaryCounterRow row, BamFile bam_file)
    
    cdef tuple _parse_column(self, JointBinaryCounterRow row, ExtendedPileupColumn pileup_column)
    
    cdef tuple _get_joint_snv_mix_features(self, JointBinaryCounterRow row)
    
    cdef tuple _get_site_features(self, ExtendedPileupColumn pileup_column)
    cdef double _get_map_qual_rms(self, ExtendedPileupColumn pileup_column)
    cdef int _get_num_zero_map_qual_reads(self, ExtendedPileupColumn pileup_column)
    cdef int _get_total_allele_count(self, ExtendedPileupColumn pileup_column, int base_qual, int map_qual)
    
    cdef tuple _get_base_features(self, char * base, ExtendedPileupColumn pileup_column)
    cdef double _get_base_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)    
    cdef double _get_map_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)    
    cdef int _get_forward_strand_count(self, char * base, ExtendedPileupColumn pileup_column, int min_base_qual)    
    cdef double _get_tail_distance_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)
    
    cdef double _get_base_matches_before_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)
    cdef double _get_base_matches_after_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)
    cdef double _get_base_matches_spanning_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)  
    
    cdef tuple _get_ratio_features(self, features)
    cdef double _get_normalised_ratio(self,
                                      int normal_depth,
                                      int tumour_depth,
                                      double normal_feature,
                                      double tumour_feature)
    
        
    
