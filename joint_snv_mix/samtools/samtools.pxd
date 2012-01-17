'''
Created on 2012-01-17

@author: Andrew Roth
'''
import os

from libc.stdint cimport int32_t, uint8_t, uint32_t, uint64_t
from libc.stdio cimport FILE
from libc.stdlib cimport free
from libc.string cimport strdup

cdef extern from "ctype.h":
    int toupper(int c)

cdef extern from * :
    ctypedef char * const_char_ptr "const char*"

cdef extern from "bam.h":
    int BAM_DEF_MASK
    
    ctypedef struct tamFile:
        pass

    ctypedef struct bamFile:
        pass        
    
    ctypedef struct bam1_core_t:
        int32_t tid 
        int32_t pos
        uint32_t bin
        uint32_t qual
        uint32_t l_qname
        uint32_t flag
        uint32_t n_cigar
        int32_t l_qseq
        int32_t mtid 
        int32_t mpos 
        int32_t isize
    
    ctypedef struct bam1_t:
        bam1_core_t core
        int l_aux
        int data_len
        int m_data
        uint8_t * data
    
    ctypedef struct bam_pileup1_t:
        bam1_t * b 
        int32_t qpos 
        int indel
        int level
        uint32_t is_del
        uint32_t is_head
        uint32_t is_tail
    
    ctypedef struct bam_header_t:
        int32_t n_targets
        char ** target_name
        uint32_t * target_len
        void * hash
        void * rg2lib
        int l_text
        char * text
    
    ctypedef struct bam_index_t:
        int32_t n
        uint64_t n_no_coor
    
    ctypedef struct bam_plbuf_t:
        pass

    ctypedef struct pair64_t:
        uint64_t u, v
            
    ctypedef struct bam_iter_t:
        int from_first
        int tid, beg, end, n_off, i, finished
        uint64_t curr_off
        pair64_t * off

    uint8_t * bam1_seq( bam1_t * b)
    uint8_t * bam1_qual( bam1_t * b)
    
    bam1_t * bam_dup1(bam1_t * src)

    bam1_t * bam_init1()
    void bam_destroy1(bam1_t *)

    bam_index_t * bam_index_load(char * f)

    void bam_index_destroy(bam_index_t * idx)
    
    int bam_parse_region(bam_header_t * header, char * str, int * ref_id, int * begin, int * end)        

    ctypedef struct bam_plp_t:
        pass
    
    ctypedef bam_pileup1_t * const_bam_pileup1_t_ptr "const bam_pileup1_t *"    
    
    ctypedef int (*bam_plp_auto_f)(void * data, bam1_t * b)
    
    bam_plp_t bam_plp_init(bam_plp_auto_f func, void * data)    
    bam_pileup1_t * bam_plp_next(bam_plp_t iter, int * _tid, int * _pos, int * _n_plp)
    bam_pileup1_t * bam_plp_auto(bam_plp_t iter, int * _tid, int * _pos, int * _n_plp)
    void bam_plp_set_mask(bam_plp_t iter, int mask)
    void bam_plp_reset(bam_plp_t iter)
    void bam_plp_destroy(bam_plp_t iter)
    
    int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t * b)    
    
    bam_iter_t bam_iter_query(bam_index_t * idx, int tid, int beg, int end)    

cdef extern from "faidx.h":
    ctypedef struct faidx_t:
        pass

    void fai_destroy(faidx_t * fai)
    
    faidx_t * fai_load(char * fn)    

    char * faidx_fetch_seq(faidx_t * fai, char * c_name, int p_beg_i, int p_end_i, int * len)

cdef extern from "sam.h":
    ctypedef struct samfile_t_un:
        tamFile tamr
        bamFile bam
        FILE * tamw
        
    ctypedef struct samfile_t:
        int type
        samfile_t_un x
        bam_header_t * header
    
    samfile_t * samopen(const_char_ptr fn, char * mode, void * aux)

ctypedef struct __iter_data:
    int tid
    int seq_len
    
    char * seq
    
    bam_iter_t iter
    
    faidx_t * fasta_file_ptr
    
    samfile_t * bam_file_ptr

#=======================================================================================================================
# Class definitions
#=======================================================================================================================
cdef class FastaFile:
    cdef char * _file_name
    cdef faidx_t * _fasta_file
   
    cdef faidx_t * get_file_pointer(self)
    cdef char * get_reference_base(self, char * reference, int position)
    cdef char * _fetch(self, char * reference, int start, int end, int * length)

cdef class BamFile:
    cdef char * _file_name
    cdef samfile_t * _bam_file
    cdef bam_index_t * _index
    
    cdef samfile_t * get_file_pointer(self)
    cdef bam_index_t * get_index(self)

cdef class AlignedRead:
    cdef  bam1_t * _delegate

cdef makeAlignedRead(bam1_t * src)

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

cdef makePileupProxy(bam_pileup1_t * plp, int tid, int pos, int n)
cdef makePileupRead(bam_pileup1_t * src)
