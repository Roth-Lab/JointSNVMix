from libc.stdlib cimport malloc, free

from csamtools cimport Samfile, Fastafile

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport PairedSampleBinomialCounterRow, get_non_ref_base, get_binary_counts
from joint_snv_mix.counters.ref_iterator cimport JointRefIterator

from joint_snv_mix.counters.quality_counter cimport QualityCounter, QualityCounterRow, QualityCounterRefIterator

from joint_snv_mix.counters.shared cimport binary_counts_struct, base_counts_struct, base_map_qualities_struct, \
                                           strcmp, toupper
   
cdef class JointBinaryQualityCounter(Counter):
    cdef QualityCounter _normal_counter
    cdef QualityCounter _tumour_counter    
    cdef Fastafile _ref_genome_fasta    
    
cdef class JointBinaryQualityCounterRow(PairedSampleBinomialCounterRow):
    cdef base_map_qualities_struct _normal_data
    cdef base_map_qualities_struct _tumour_data

cdef class JointBinaryQualityCounterIterator(JointRefIterator): 
    cdef Fastafile _ref_genome_fasta
