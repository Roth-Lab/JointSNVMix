from libc.stdlib cimport free

from csamtools cimport Samfile, Fastafile

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport PairedSampleCounterRow, get_non_ref_base, get_binary_counts
from joint_snv_mix.counters.ref_iterator cimport JointRefIterator
from joint_snv_mix.counters.base_counter cimport BaseCounter, BaseCounterRow, BaseCounterRefIterator
from joint_snv_mix.counters.shared cimport binary_counts_struct, base_counts_struct, strcmp, counts_struct
   
cdef class JointBinaryBaseCounter(Counter):
    cdef BaseCounter _normal_counter
    cdef BaseCounter _tumour_counter    
    cdef Fastafile _ref_genome_fasta    
    
cdef class JointBinaryCounterRow(PaireSampleBinomialCounterRow):
    cdef binary_counts_struct _normal_counts
    cdef binary_counts_struct _tumour_counts

cdef class JointBinaryBaseCounterIterator(JointRefIterator):
    cdef Fastafile _ref_genome_fasta
