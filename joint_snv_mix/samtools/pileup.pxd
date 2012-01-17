'''
Created on 2012-01-17

@author: Andrew Roth
'''
from libc.stdint cimport int32_t, uint8_t, uint32_t, uint64_t
from libc.stdlib cimport free

from joint_snv_mix.samtools.samtools_clib import bam_destroy1, bam_iter_read, bam_plp_reset, bam_plp_init, \
                                                 bam_plp_destroy, bam_plp_auto, bam_plp_set_mask, BAM_DEF_MASK, \
                                                 bam_dup1, bam_iter_t, bam1_t, samfile_t, const_bam_pileup1_t_ptr, \
                                                 bam_plp_t, bam_pileup1_t

cdef class AlignedRead:
    cdef  bam1_t * _delegate

cdef class IteratorRowRegion(object):
    cdef int _return_value    
    
    cdef bam_iter_t _iter
    
    cdef bam1_t * _bam_struct
    cdef samfile_t * _bam_file_ptr
    
    cdef BamFile _bam_file

    cdef bam1_t * getCurrent(self)
    cdef cnext(self)

cdef class IteratorColumnRegion:
    cdef int _tid
    cdef int _pos
    cdef int _n_plp
    cdef int _mask
    
    cdef BamFile _bam_file
#    cdef FastaFile _fasta_file
    
    cdef const_bam_pileup1_t_ptr _plp    
    cdef bam_plp_t _pileup_iter    
    cdef __iter_data _iter_data 
    
    cdef IteratorRowRegion _iter

    cdef cnext(self)
    cdef _setup_iterator_data(self, int tid, int start, int end)
    
cdef class PileupProxy:
    cdef bam_pileup1_t * _plp
    cdef int _tid
    cdef int _pos
    cdef int _n_pu
   
cdef class PileupRead:
    cdef int _indel
    cdef int _level
    
    cdef int32_t  _qpos
   
    cdef uint32_t _is_del
    cdef uint32_t _is_head
    cdef uint32_t _is_tail
    
    cdef AlignedRead _alignment

cdef makeAlignedRead(bam1_t * src)
cdef makePileupProxy(bam_pileup1_t * plp, int tid, int pos, int n)
cdef makePileupRead(bam_pileup1_t * src)
