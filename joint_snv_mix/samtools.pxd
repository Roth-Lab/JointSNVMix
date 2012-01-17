'''
Created on 2012-01-17

@author: Andrew Roth
'''
cdef class Fastafile:
    cdef char * _file_name
    cdef faidx_t * _fasta_file
   
    cdef faidx_t * get_file_pointer(self)
    cdef char * get_reference_base(self, char * reference, int position)
    cdef char * _fetch(self, char * reference, int start, int end, int * length)

cdef class BamFile:
    cdef char * _file_name
    cdef samfile_t * _bam_file
    cdef bam_index_t * _index

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
    cdef int cnext(self)

cdef class IteratorColumn:
    cdef int _tid
    cdef int _pos
    cdef int _n_plp
    cdef int _mask
    
    cdef BamFile _bam_file
    cdef FastaFile _fasta_file
    
    cdef const_bam_pileup1_t_ptr _plp    
    cdef bam_plp_t _pileup_iter    
    cdef __iterdata _iter_data 

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
