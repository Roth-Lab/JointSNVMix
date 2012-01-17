'''
Created on 2012-01-17

@author: Andrew Roth
'''
cdef class AlignedRead:
    def __dealloc__(self):
        bam_destroy1(self._delegate)

cdef class IteratorRowRegion(object):
    def __cinit__(self, BamFile bam_file, int tid, int start, int end):        
        # Makes sure that samfile stays alive as long as the iterator is alive.
        self._bam_file = bam_file

        self._bam_file_ptr = self._bam_file.get_file_pointer()

        self._return_value = 0

        self._iter = bam_iter_query(self._bam_file.get_index(), tid, start, end)

        self._bam_struct = bam_init1()
        
    def __dealloc__(self):
        bam_destroy1(self._bam_struct)        

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent(self):
        return self._bam_struct

    cdef cnext(self):
        self._return_value = bam_iter_read(self._bam_file_ptr.x.bam,
                                           self._iter,
                                           self._bam_struct)
        
    def __next__(self): 
        self.cnext()

        if self.retval < 0: 
            raise StopIteration
        
        return makeAlignedRead(self._bam_struct)

cdef class IteratorColumnRegion:
    def __cinit__(self, BamFile bam_file, int tid=0, int start=0, int end=max_pos):        
        self._bam_file = bam_file
        
#        self._fasta_file = fasta_file
                                
        self._mask = BAM_DEF_MASK

        self._tid = 0
        self._pos = 0
        self._n_plp = 0
        self._plp = NULL
                
        self._setup_iterator_data(self._tid, start, end)

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators that have
        # not been fully consumed
        if self._pileup_iter != < bam_plp_t > NULL:
            bam_plp_reset(self._pileup_iter)
            
            bam_plp_destroy(self._pileup_iter)
            
            self._pileup_iter = < bam_plp_t > NULL
    
        if self._iter_data.seq != NULL: 
            free(self._iter_data.seq)
            
            self._iter_data.seq = NULL
        
    def __iter__(self):
        return self

    def __next__(self): 
        while True:
            self.cnext()
            
            if self._n_plp < 0:
                raise Exception("IteratorColumn : error during iteration.")
        
            if self._plp == NULL:
                raise StopIteration
            
            return makePileupProxy(< bam_pileup1_t *> self._plp,
                                     self._tid,
                                     self._pos,
                                     self._n_plp)    
    
    cdef cnext(self):
        '''
        Perform next iteration.
    
        This method is analogous to the samtools bam_plp_auto method.
        It has been re-implemented to permit for filtering.
        '''
        self._plp = bam_plp_auto(self._pileup_iter,
                                 & self._tid,
                                 & self._pos,
                                 & self._n_plp)
        
    cdef _setup_iterator_data(self, int tid, int start, int end):
        '''
        Setup the iterator structure.
        '''    
        self._iter = IteratorRowRegion(self._bam_file, tid, start, end)
        
        self._iter_data.bam_file_ptr = self._bam_file.get_file_pointer()
        self._iter_data.fasta_file_ptr = NULL
        self._iter_data.iter = self._iter._iter
        self._iter_data.seq = NULL
        self._iter_data.tid = -1
            
        self._pileup_iter = bam_plp_init(& __advance_all, & self._iter_data)
    
        bam_plp_set_mask(self._pileup_iter, self._mask)

#=======================================================================================================================
# Pileup objects
#=======================================================================================================================
cdef class PileupProxy:
    pass

cdef class PileupRead:
    pass

#=======================================================================================================================
# Factory methods
#=======================================================================================================================
cdef makeAlignedRead(bam1_t * src):
    cdef AlignedRead dest = AlignedRead.__new__(AlignedRead)
    
    dest._delegate = bam_dup1(src)
    
    return dest
 
cdef makePileupProxy(bam_pileup1_t * plp, int tid, int pos, int n):
    cdef PileupProxy dest = PileupProxy.__new__(PileupProxy)
    
    dest._plp = plp
    dest._tid = tid
    dest._pos = pos
    dest._n_pu = n
    
    return dest

cdef makePileupRead(bam_pileup1_t * src):
    cdef PileupRead dest = PileupRead.__new__(PileupRead)
    
    dest._alignment = makeAlignedRead(src.b)
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    
    return dest
        
#=======================================================================================================================
# Utility functions and structs
#=======================================================================================================================
cdef int __advance_all(void * data, bam1_t * b):
    '''
    Advance without any read filtering.
    '''
    cdef __iter_data * d
    
    d = < __iter_data *> data
    
    return bam_iter_read(d.bam_file_ptr.x.bam, d.iter, b)
