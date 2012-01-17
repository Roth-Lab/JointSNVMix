'''
Created on 2011-06-21

@author: Andrew Roth
'''
cdef class BaseCounter(Counter):
    '''
    Class for counting all four bases at each position in a bam file.
    '''
    def __init__(self, BamFile bam_file, int min_base_qual=10, int min_map_qual=10):
        self._bam_file = bam_file
        
        self._min_base_qual = min_base_qual
        self._min_map_qual = min_map_qual
        
        self._refs = self._bam_file.references
    
    def iter_ref(self, ref):
        '''
        Returns a BaseCounterIter iterator for given ref.
        '''
        iter = BaseCounterRefIterator(ref, self._bam_file.pileup(ref), self._min_base_qual, self._min_map_qual)
        
        return iter

cdef class BaseCounterRefIterator(RefIterator):
    '''
    Iterator class for iterating over a reference in a bam file and returning BaseCounterRow objects.
    '''
    def __init__(self, char * ref, IteratorColumnRegion pileup_iter, int min_base_qual, int min_map_qual):
        self._ref = ref
                
        self._min_base_qual = min_base_qual
        self._min_map_qual = min_map_qual
        
        self._ref_iter = CRefIterator(ref, pileup_iter)                                
        self._position = self._ref_iter._position
  
    cdef advance_position(self):
        self._ref_iter.advance_position()
        self._position = self._ref_iter._position
    
    cdef parse_current_position(self):
        cdef column_struct column
        cdef counts_struct counts
                
        self._ref_iter.parse_current_position()
        
        column = self._ref_iter._current_column        
        
        counts = self._parse_column(column)
    
        self._current_row = makeBaseCounterRow(self._ref, self._position, counts)    
            
    cdef counts_struct _parse_column(self, column_struct column):            
        cdef int i
        cdef char base[2]
        cdef counts_struct counts
        
        counts.A = 0
        counts.C = 0
        counts.G = 0
        counts.T = 0
        
        for i in range(column.depth):                     
            if column.map_quals[i] < self._min_map_qual:
                continue            
                                 
            if column.base_quals[i] < self._min_base_qual:
                continue
            
            base[0] = column.bases[i]
            base[1] = '\0'
            
            if strcmp(base, "A") == 0:
                counts.A += 1        
            elif strcmp(base, "C") == 0:
                counts.C += 1
            elif strcmp(base, "G") == 0:
                counts.G += 1
            elif strcmp(base, "T") == 0:
                counts.T += 1
        
        return counts

cdef class BaseCounterRow(SingleSampleCounterRow):
    '''
    Class for storing four nucleotide count data from Bam file position.
    '''
    def __str__(self):
        row = [self.ref, str(self.position)]
        row.extend([str(x) for x in self.counts])
        
        return "\t".join(row)
        
    property counts:
        '''
        Return the counts for the four bases as tuple in the order A,C,G,T.
        '''
        def __get__(self):
            return (
                    self._counts.A,
                    self._counts.C,
                    self._counts.G,
                    self._counts.T
                    )

    cdef int get_counts(self, char * base):
        cdef int result
    
        if strcmp(base, "A") == 0:
            result = self._counts.A        
        elif strcmp(base, "C") == 0:
            result = self._counts.C
        elif strcmp(base, "G") == 0:
            result = self._counts.G
        elif strcmp(base, "T") == 0:
            result = self._counts.T
        else:
            return 0
        
        return result
'''
C level constructor for BaseCounterRow object.
'''
cdef class BaseCounterRow
cdef BaseCounterRow makeBaseCounterRow(char * ref, int position, counts_struct counts):
    cdef BaseCounterRow row = BaseCounterRow.__new__(BaseCounterRow)
    
    row._ref = ref
    row._position = position
    row._counts = counts     
    row._depth = counts.A + counts.C + counts.G + counts.T
    
    return row
