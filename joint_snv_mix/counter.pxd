'''
Classes for creating iterators for count data over a pair of genomes. 

Created on 2012-01-18

@author: Andrew Roth
'''
cdef class JointBinaryCounter(Counter):
        cdef char * _type
        
        cdef BamFile _normal_bam
        cdef BamFile _tumour_bam
        
        cdef FastaFile _ref_genome
        
        cdef tuple _refs 
        
cdef class JointBinaryCounterIterator(RefIterator):
    cdef char * _type
    cdef char * _ref

    cdef int _min_base_qual
    cdef int _min_map_qual
    
    cdef PileupIterator _normal_iter
    cdef PileupIterator _tumour_iter
    
    cdef FastaFile _ref_genome
    
    cdef _int pos
        
    cdef advance_position(self)
    cdef parse_current_position(self)
    
    cdef _make_counter_row(self, normal_column, tumour_column)
    cdef _make_count_data(self, ref_base, var_base, normal_column, tumour_column)
    cdef _make_quality_data(self, ref_base, normal_column, tumour_column)

cdef class JointBinaryCounterRow(object):
    cdef char * _ref
    cdef int _pos
    
    cdef char * _ref_base
    cdef char * _var_base
    
    cdef JointBinaryData _data
    
cdef class JointBinaryBaseCounterRow(JointBinaryCounterRow):
    pass

cdef class JointBinaryQualityCounterRow(JointBinaryCounterRow):
    pass

#=======================================================================================================================
# Data object 
#=======================================================================================================================
cdef class JointBinaryData(object):
    pass

cdef class JointBinaryCountData(object):
    cdef int _a_N
    cdef int _b_N
    
    cdef int _a_T
    cdef int _b_T   

cdef class JointBinaryQualityData(object):
    cdef int _normal_depth
    cdef int _tumour_depth
    
    cdef double * _q_N
    cdef double * _r_N
    
    cdef double * _q_T
    cdef double * _r_T