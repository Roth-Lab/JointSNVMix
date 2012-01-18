'''
Created on 2012-01-17

@author: Andrew Roth
'''
# Global constants
cdef int max_pos = 2 << 29
cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"

#cdef class AlignedRead:
#    def __dealloc__(self):
#        bam_destroy1(self._delegate)

cdef class PileupRegionIterator:
    def __cinit__(self, BamFile bam_file, int tid=0, int start=0, int stop=max_pos):        
        self._bam_file = bam_file
                            
        self._mask = BAM_DEF_MASK

        self._tid = tid
        self._pos = 0
        self._n_plp = 0
        self._plp = NULL
                
        self._setup_iterator_data(self._tid, start, stop)

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
            
            return makePileupColumn(< bam_pileup1_t *> self._plp,
                                    self._tid,
                                    self._pos,
                                    self._n_plp)
#            return self._tid, self._pos, self._n_plp
    
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
        
    cdef _setup_iterator_data(self, int tid, int start, int stop):
        '''
        Setup the iterator structure.
        '''
        cdef bam_iter_t iter
        
        iter = bam_iter_query(self._bam_file.get_index(), tid, start, stop)
        
        self._iter_data.bam_file_ptr = self._bam_file.get_file_pointer()
        self._iter_data.iter = iter
        self._iter_data.seq = NULL
        self._iter_data.tid = -1
            
        self._pileup_iter = bam_plp_init(& advance_all, & self._iter_data)
    
        bam_plp_set_mask(self._pileup_iter, self._mask)
        
        print self._iter_data.tid

#cdef class PileupProxy:
#    pass

cdef class PileupColumn:
    def __dealloc__(self):
        free(self._bases)
        free(self._base_quals)
        free(self._map_quals)
    
    def __str__(self):
        return self._tid, self._pos

#cdef class PileupRead:
#    pass

#=======================================================================================================================
# Factory methods
#=======================================================================================================================
#cdef makeAlignedRead(bam1_t * src):
#    cdef AlignedRead dest = AlignedRead.__new__(AlignedRead)
#    
#    dest._delegate = bam_dup1(src)
#    
#    return dest
 
#cdef makePileupProxy(bam_pileup1_t * plp, int tid, int pos, int n):
#    cdef PileupProxy dest = PileupProxy.__new__(PileupProxy)
#    
#    dest._num_reads = n
#    dest._pileup_ptr = plp
#    
#    dest._tid = tid
#    dest._pos = pos
#
#    return dest

#cdef makePileupRead(bam_pileup1_t * src):
#    cdef PileupRead dest = PileupRead.__new__(PileupRead)
#    
#    dest._alignment = makeAlignedRead(src.b)
#    
#    dest._query_pos = src.qpos
#    
#    dest._indel_size = src.indel
#
#    dest._level_in_viewer = src.level
#
#    dest._is_del = src.is_del
#    dest._is_head = src.is_head
#    dest._is_tail = src.is_tail
#    
#    return dest

cdef makePileupColumn(bam_pileup1_t * plp, int tid, int pos, int n):
    cdef int i, depth, index, qpos
    cdef bam1_t * alignment
    cdef bam_pileup1_t * pileup

    cdef PileupColumn column = PileupColumn.__new__(PileupColumn)
                    
    column._tid = tid
    column._pos = pos
    
    # Find out how many non deletion reads there are.
    depth = get_covered_depth(plp, n)    
    
    column._depth = depth
     
    column._bases = < char *> malloc(sizeof(char) * depth)
    column._base_quals = < int *> malloc(sizeof(int) * depth)
    column._map_quals = < int *> malloc(sizeof(int) * depth)  
        
    index = 0

    for i in range(n):
        pileup = & plp[i]
        
        if pileup.is_del:
            continue
        
        qpos = pileup.qpos
        
        alignment = bam_dup1(pileup.b)
        
        column._bases[index] = get_base(alignment, qpos)            

        column._base_quals[index] = get_qual(alignment, qpos)
        
        column._map_quals[index] = alignment.core.qual
        
        bam_destroy1(alignment)
        
        index += 1
    
    return column
        
#=======================================================================================================================
# Utility functions and structs
#=======================================================================================================================
cdef int advance_all(void * data, bam1_t * b):
    '''
    Advance without any read filtering.
    '''
    cdef iter_data_t * d
    
    d = < iter_data_t *> data
    
    return bam_iter_read(d.bam_file_ptr.x.bam, d.iter, b)

cdef char get_base(bam1_t * src, int pos):
    cdef uint8_t * p
    cdef char base

    if src.core.l_qseq == 0: 
        return None
    
    if not src.core.l_qseq:
        return None

    seq = bam1_seq(src)
    
    base = bam_nt16_rev_table[seq[pos / 2] >> 4 * (1 - pos % 2) & 0xf]
    
    return base

cdef int get_qual(bam1_t * src, int pos):    
    cdef uint8_t * p

    p = bam1_qual(src)
    if p[0] == 0xff:
        return None

    return p[pos]

cdef int get_covered_depth(bam_pileup1_t * plp, int n):
    '''
    Gets the number of reads without deletions in a proxy object.
    '''
    cdef int i, depth
    cdef bam_pileup1_t * pileup
    
    depth = 0
    
    for i in range(n):
        pileup = & plp[i]
        
        if pileup.is_del:
            continue
        
        depth += 1
    
    return depth
