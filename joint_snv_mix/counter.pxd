'''
Classes for creating iterators for count data over a pair of genomes. 

Created on 2012-01-18

@author: Andrew Roth
'''
from joint_snv_mix.counter_row cimport RowFactory
from joint_snv_mix.ref_iterator cimport RefIterator

from joint_snv_mix.samtools.bam cimport BamFile
from joint_snv_mix.samtools.fasta cimport FastaFile
from joint_snv_mix.samtools.pileup cimport PileupColumn, PileupIterator


cdef class JointBinaryCounter(object):
    cdef bint _qualities
    
    cdef int _min_base_qual
    cdef int _min_map_qual
    
    cdef BamFile _normal_bam
    cdef BamFile _tumour_bam
    
    cdef FastaFile _ref_genome
    
    cdef tuple _refs 
        
cdef class JointBinaryCounterIterator(RefIterator):    
    cdef PileupIterator _normal_iter
    cdef PileupIterator _tumour_iter
    
    cdef RowFactory _row_factory
    
    cdef jump_to_position(self, int position)