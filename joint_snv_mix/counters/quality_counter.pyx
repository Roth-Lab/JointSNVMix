'''
Use the convention every position at C-level is 0-based. If the position is
accessed by property mechanism from python layer it should be one based.

Created on 2011-07-28

@author: Andrew Roth
'''
DEF NUM_QUAL_VAL = 256

cdef class QualityCounter(Counter):
    '''
    Class for counting all four bases at each position in a bam file.
    '''
    def __init__(self, Samfile bam_file):
        self._bam_file = bam_file
        
        self._refs = self._bam_file.references
    
    def iter_ref(self, ref):
        '''
        Returns a BaseCounterIter iterator for given ref.
        '''
        iter = QualityCounterRefIterator(ref, self._bam_file.pileup(ref))
        
        return iter

cdef class QualityCounterRefIterator(CounterRefIterator):
    '''
    Iterator class for iterating over a reference in a bam file and returning BaseCounterRow objects.
    '''
    def __init__(self, char * ref, IteratorColumnRegion pileup_iter):
        self._ref = ref
        
        self._ref_iter = CRefIterator(ref, pileup_iter)                                
        self._position = self._ref_iter._position
  
    cdef cnext(self):
        self.advance_position()
        self.parse_current_position()
    
    cdef advance_position(self):
        self._ref_iter.advance_position()
        self._position = self._ref_iter._position
    
    cdef parse_current_position(self):
        cdef column_struct column
                
        self._ref_iter.parse_current_position()
        
        column = self._ref_iter._current_column        
    
        self._current_row = makeQualityCounterRow(column)  
    
    
cdef class QualityCounterRow(CounterRow):
    '''
    Class for storing count data from Bam file position.
    '''    
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")
    
    def __dealloc__(self):
        free(self._bases)
        free(self._base_quals)
        free(self._map_quals)
    
    property counts:
        '''
        Return the counts for the four bases as tuple in the order A,C,G,T.
        '''
        def __get__(self):
            cdef char base[2]
            cdef list counts = []
            
            base[1] = '\0'
            
            for i in range(self._depth):
                base[0] = self._bases[i]                
                
                counts.append((
                               base,
                               self._base_quals[i],
                               self._map_quals[i]
                               ))
            
            
            return counts
                
    property depth:
        '''
        Depth of all counts.
        '''
        def __get__(self):            
            return self._depth

    cdef base_counts_struct get_base_counts(self, char * base):
        cdef int i
        cdef base_counts_struct base_counts
        
        base_counts.base = base
        base_counts.counts = 0
        
        for i in range(self._depth):
            if base[0] == self._bases[i]:
                base_counts.counts += 1
        
        return base_counts
    
    cdef int get_counts(self, char * base):
        cdef int i, counts
        
        counts = 0
        
        for i in range(self._depth):
            if base[0] == self._bases[i]:
                counts += 1
        
        return counts

'''
C level constructor for BaseCounterRow object.
'''
cdef class QualityCounterRow
cdef QualityCounterRow makeQualityCounterRow(column_struct column):
    cdef int i

    cdef QualityCounterRow row = QualityCounterRow.__new__(QualityCounterRow)
    
    row._ref = column.ref
    row._position = column.position     
    row._depth = column.depth
    
    row._bases = < char *> malloc(column.depth * sizeof(char))
    
    row._base_quals = < int *> malloc(column.depth * sizeof(int))
    
    row._map_quals = < int *> malloc(column.depth * sizeof(int))
    
    for i in range(column.depth):
        row._bases[i] = column.bases[i]
        row._base_quals[i] = column.base_quals[i]
        row._map_quals[i] = column.map_quals[i]     
    
    return row
                         
#===============================================================================
# Modified pysam code
#===============================================================================
cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"

cdef char * get_base(bam1_t * src, int pos):
    cdef uint8_t * p
    cdef char base
    cdef char[2] base_str

    if src.core.l_qseq == 0: 
        return None
    
    if not src.core.l_qseq:
        return None

    seq = bam1_seq(src)
    
    base = bam_nt16_rev_table[seq[pos / 2] >> 4 * (1 - pos % 2) & 0xf]
    
    base_str[0] = base
    base_str[1] = '\0'
    
    return base_str

cdef int get_qual(bam1_t * src, int pos):    
    cdef uint8_t * p

    p = bam1_qual(src)
    if p[0] == 0xff:
        return None

    return p[pos]
