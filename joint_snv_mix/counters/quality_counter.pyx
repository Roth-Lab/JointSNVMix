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

cdef class QualityCounterRefIterator(RefIterator):
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
    
    
cdef class QualityCounterRow(SingleSampleCounterRow):
    def __dealloc__(self):
        free(self._bases)
        free(self._base_quals)
        free(self._map_quals)
    
    def __str__(self):
        row = [self.ref, str(self.position)]
        
        return "\t".join(row)
    
    property bases:        
        def __get__(self):
            cdef char base[2]
            base[1] = '\0'
            
            bases = []

            for i in range(self._depth):
                base[0] = self._bases[i]                                
                bases.append(base)            
            
            return bases
    
    property map_quals:        
        def __get__(self):
            map_quals = []
            
            for i in range(self._depth):
                map_quals.append(self._map_quals[i])                                
            
            return map_quals
    
    property base_quals:        
        def __get__(self):
            base_quals = []
            
            for i in range(self._depth):
                base_quals.append(self._base_quals[i])                                
            
            return base_quals
    
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
