'''
Created on 2012-01-26

@author: Andrew Roth
'''
from libc.math cimport pow

from joint_snv_mix.samtools.bam cimport BamFile
from joint_snv_mix.samtools.pileup cimport ExtendedPileupColumn

from joint_snv_mix.counter_row cimport JointBinaryCounterRow

cdef class FeatureExtractor(object):
    cdef BamFile _normal_bam
    cdef BamFile _tumour_bam
   
    cdef list _get_genome_features(self, JointBinaryCounterRow row, BamFile bam_file)
    
    cdef list _parse_column(self, JointBinaryCounterRow row, ExtendedPileupColumn pileup_column)
    
    cdef list _get_base_features(self, char * base, ExtendedPileupColumn pileup_column)
    
    cdef double _get_base_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)
    
    cdef double _get_map_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)
    
    cdef int _get_forward_strand_count(self, char * base, ExtendedPileupColumn pileup_column)
    
    cdef double _get_tail_distance_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent)         
    
            
        
    
