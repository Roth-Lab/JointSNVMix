cdef class CRefIterator(object):
    '''
    Designed for C level access only. Main attribute is current_column which is column_struct object.
    
    Iterator class for iterating over a reference in a bam file data for each position.
    '''
    def __init__(self, char * ref, IteratorColumnRegion pileup_iter):
        self._ref = ref
        self._position = -1
        
        self._pileup_iter = pileup_iter
        
        self._current_column.bases = NULL
    
    def __dealloc__(self):
        if self._current_column.bases is not NULL:
            destroy_column_struct(self._current_column)
  
    cdef cnext(self):
        '''
        C-level iterator function. Moves current_column to next position.
        '''
        self.advance_position()
        self.parse_current_position()
        
    cdef advance_position(self):
        '''
        Advance pileup iter without parsing the column. Only self.position attribute of object should be examined,
        self.current_column is in an undefined state unless self.parse_current_position() is called.
        '''
        self._current_pileup_column = self._pileup_iter.next()        
        self._position = self._current_pileup_column.pos
    
    cdef parse_current_position(self):
        '''
        Parses current value stored in self.current_pileup_column. 
        '''
        if self._current_column.bases is not NULL:
            destroy_column_struct(self._current_column)
        
        self._current_column = self._parse_pileup_column(self._current_pileup_column)
                
    cdef column_struct _parse_pileup_column(self, PileupProxy pileup_column):
        cdef int i, depth, index, qpos
        cdef bam1_t * alignment
        cdef bam_pileup1_t * pileup
        cdef column_struct column
                
        depth = self._get_depth(pileup_column)        
        column = make_column_struct(self._ref, pileup_column.pos, depth)
        
        index = 0

        for i in range(pileup_column.n_pu):
            pileup = & pileup_column.plp[i]
            
            if pileup.is_del:
                continue
            
            qpos = pileup.qpos
            
            alignment = bam_dup1(pileup.b)
            
            column.bases[index] = get_base(alignment, qpos)            

            column.base_quals[index] = get_qual(alignment, qpos)
            
            column.map_quals[index] = alignment.core.qual
            
            bam_destroy1(alignment)
            
            index += 1
            
        return column
    
    cdef int _get_depth(self, PileupProxy pileup_column):
        cdef int i, depth
        
        depth = 0
        
        for i in range(pileup_column.n_pu):
            pileup = & pileup_column.plp[i]
            
            if pileup.is_del:
                continue
            
            depth += 1
        
        return depth
    
cdef column_struct make_column_struct(char * ref, int position, int depth):
    cdef column_struct column
    
    column.ref = ref
    
    column.position = position
    
    column.depth = depth
    
    column.bases = < char *> malloc(depth * sizeof(char))
    
    column.base_quals = < int *> malloc(depth * sizeof(int))
    
    column.map_quals = < int *> malloc(depth * sizeof(int))
    
    return column

cdef void destroy_column_struct(column_struct column):
    free(column.bases)
    free(column.base_quals)
    free(column.map_quals)

#=======================================================================================================================
# Modified pysam code
#=======================================================================================================================
cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"

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
