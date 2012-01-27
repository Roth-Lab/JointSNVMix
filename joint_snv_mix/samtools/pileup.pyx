'''
Created on 2012-01-17

@author: Andrew Roth
'''
# Global constants
cdef int max_pos = 2 << 29
cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"

cdef class PileupIterator:
    def __cinit__(self, BamFile bam_file, int tid=0, int start=0, int stop=max_pos):        
        self._bam_file = bam_file
                            
        self._mask = BAM_DEF_MASK

        self._tid = tid
        self._pos = 0
        self._n_plp = 0
        self._plp = NULL
                
        self._setup_iterator_data(self._tid, start, stop)

    def __dealloc__(self):
        self._destroy_iterator_data()
        
    def __iter__(self):
        return self

    def __next__(self): 
        self.cnext()

        return self._column
    
    cdef cnext(self):
        '''
        C-level iterator function. Moves current_column to next position.
        '''
        self.advance_position()
        self.parse_current_position()
    
    cdef advance_position(self):
        '''
        This function moves the underlying pileup pointer along.
        '''
        self._plp = bam_plp_auto(self._pileup_iter,
                                 & self._tid,
                                 & self._pos,
                                 & self._n_plp)
        
        if self._plp == NULL:
            raise StopIteration
    
    cdef parse_current_position(self):
        '''
        This function parses the current position to create a pileup column.
        '''
        if self._n_plp < 0:
            raise Exception("IteratorColumn : error during iteration.")        

        self._column = makePileupColumn(< bam_pileup1_t *> self._plp,
                                        self._tid,
                                        self._pos,
                                        self._n_plp)
    
    cpdef jump_to_position(self, int position):
        '''
        Move iterator to position. Should probably only be used for big jumps since there is some overhead from
        destroying and constructing underlying c data structues.
        '''        
        self._destroy_iterator_data()
        self._setup_iterator_data(self._tid, position-1, max_pos)
        self.advance_position()
        
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
    
    cdef _destroy_iterator_data(self):
        # reset in order to avoid memory leak messages for iterators that have
        # not been fully consumed
        if self._pileup_iter != < bam_plp_t > NULL:
            bam_plp_reset(self._pileup_iter)
            
            bam_plp_destroy(self._pileup_iter)
            
            self._pileup_iter = < bam_plp_t > NULL
    
        if self._iter_data.seq != NULL: 
            free(self._iter_data.seq)            
            self._iter_data.seq = NULL 
        
        if self._iter_data.iter != < bam_iter_t > NULL:
            bam_iter_destroy(self._iter_data.iter)
            self._iter_data.iter = < bam_iter_t > NULL

cdef class PileupColumn:
    def __dealloc__(self):
        free(self._bases)
        free(self._base_quals)
        free(self._map_quals)
    
    def __str__(self):
        return "{0}, {1}, {2}".format(self._tid, self._pos, self._depth)
    
    property position:
        def __get__(self):
            return self._pos + 1
    
    property base_quals:
        def __get__(self):
            return [x for x in self._base_quals[:self._depth]]

    property map_quals:
        def __get__(self):
            return [x for x in self._map_quals[:self._depth]]
    
    cdef int get_depth(self):
        return self._depth
    
    cdef int get_nucleotide_count(self, char * base, int min_base_qual, int min_map_qual):
        '''
        Retrieve the number of nucleotides found at this position. Counts will be initialised to 0 in the method.
        
        Arguments:
        base - a char * containing the query base.
        min_base_qual - nucleotides will only be counted if their base quality exceeds this value.
        min_map_qual - nucleotides will only be counted if their mapping quality exceeds this value.
        '''
        cdef int i
        cdef int count
        
        count = 0
        
        for i in range(self._depth):            
            if self._base_quals[i] < min_base_qual or self._map_quals[i] < min_map_qual:
                continue            
            
            if base[0] == self._bases[i]:
                count += 1
        
        return count
    
#=======================================================================================================================
# Factory methods
#=======================================================================================================================
cdef makePileupColumn(bam_pileup1_t * plp, int tid, int pos, int n_plp):
    cdef int i, depth, index, qpos
    cdef bam1_t * alignment
    cdef bam_pileup1_t * pileup

    cdef PileupColumn column = PileupColumn.__new__(PileupColumn)
                    
    column._tid = tid
    
    column._pos = pos
    
    column._depth = get_covered_depth(plp, n_plp)
     
    column._bases = < char *> malloc(sizeof(char) * column._depth)
    column._base_quals = < int *> malloc(sizeof(int) * column._depth)
    column._map_quals = < int *> malloc(sizeof(int) * column._depth)  
        
    index = 0

    for i in range(n_plp):
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
