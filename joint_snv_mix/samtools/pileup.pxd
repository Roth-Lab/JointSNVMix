'''
Created on 2012-01-17

@author: Andrew Roth
'''
from libc.stdint cimport int32_t, uint8_t, uint32_t, uint64_t
from libc.stdlib cimport free, malloc
from libc.string cimport strcmp, strdup

from joint_snv_mix.samtools.samtools_clib cimport bam_destroy1, bam_iter_read, bam_plp_reset, bam_plp_init, \
                                                 bam_plp_destroy, bam_plp_auto, bam_plp_set_mask, BAM_DEF_MASK, \
                                                 bam_dup1, bam_iter_t, bam1_t, samfile_t, const_bam_pileup1_t_ptr, \
                                                 bam_plp_t, bam_pileup1_t, bam_iter_query, bam1_seq, bam1_qual, bam_iter_destroy 

from joint_snv_mix.samtools.bam cimport BamFile
                                                 
ctypedef struct iter_data_t:
    int tid
    samfile_t * bam_file_ptr
    bam_iter_t iter
    int seq_len
    char * seq

cdef class PileupIterator:
    cdef int _tid
    cdef int _pos
    cdef int _n_plp
    cdef int _mask
    
    cdef BamFile _bam_file
    
    cdef const_bam_pileup1_t_ptr _plp    
    cdef bam_plp_t _pileup_iter    
    cdef iter_data_t _iter_data 
    
    cdef PileupColumn _column

    cdef cnext(self)
    
    cdef advance_position(self)
    cdef parse_current_position(self)
    
    cpdef jump_to_position(self, int position)
    
    cdef _setup_iterator_data(self, int tid, int start, int end)
    
    cdef _destroy_iterator_data(self)
        
cdef class PileupColumn:
    cdef int _tid
    cdef int _pos
    cdef int _depth
    
    cdef char * _bases
    
    cdef int * _base_quals
    cdef int * _map_quals
    
    cdef int get_depth(self)
    
    cdef int get_nucleotide_count(self, char * base, int min_base_qual, int min_map_qual)
