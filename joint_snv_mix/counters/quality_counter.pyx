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
        self._pileup_iter = pileup_iter
        
        self._position = pileup_iter.pos
        
        self._init_qual_map()
    
    cdef _init_qual_map(self):
        cdef double base, exp
        
        for qual in range(NUM_QUAL_VAL):
            exp = -1 * (< double > qual) / 10
            base = 10
            
            self._qual_map[qual] = 1 - pow(base, exp)
  
    cdef cnext(self):
        cdef PileupProxy pileup_column
        
        pileup_column = self._pileup_iter.next()
        
        self._position = pileup_column.pos
        
        self._current_row = self._parse_pileup_column(pileup_column)
            
    cdef QualityCounterRow _parse_pileup_column(self, PileupProxy pileup_column):
        cdef int i, qpos, map_qual, base_qual, num_reads
        cdef char * base        
        cdef bam1_t * alignment
        cdef bam_pileup1_t * pileup        
        cdef char * bases
        cdef double * base_quals, * map_quals
        
        num_reads = 0
        
        for i in range(pileup_column.n_pu):
            pileup = & pileup_column.plp[i]
            
            if pileup.is_del:
                continue
            
            num_reads += 1
                    
        bases = < char *> malloc(num_reads * sizeof(char))
        base_quals = < double *> malloc(num_reads * sizeof(double))
        map_quals = < double *> malloc(num_reads * sizeof(double))
        
        for i in range(num_reads):
            pileup = & pileup_column.plp[i]
            
            if pileup.is_del:
                continue
            
            qpos = pileup.qpos
            
            alignment = bam_dup1(pileup.b) 
            
            map_qual = alignment.core.qual
            map_quals[i] = self._qual_map[map_qual]
                        
            base_qual = get_qual(alignment, qpos)
            base_quals[i] = self._qual_map[base_qual]
                       
            base = get_base(alignment, qpos)
            bases[i] = base[0]

            bam_destroy1(alignment)
        
        return makeQualityCounterRow(self._ref, self._position, num_reads, bases, base_quals, map_quals)
    
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
            
            for i in range(self._num_reads):
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
            return self._num_reads

    cdef base_counts_struct get_base_counts(self, char * base):
        cdef int i
        cdef base_counts_struct base_counts
        
        base_counts.base = base
        base_counts.counts = 0
        
        for i in range(self._num_reads):
            if base[0] == self._bases[i]:
                base_counts.counts += 1
        
        return base_counts
    
    cdef int get_counts(self, char * base):
        cdef int i, counts
        
        counts = 0
        
        for i in range(self._num_reads):
            if base[0] == self._bases[i]:
                counts += 1
        
        return counts

'''
C level constructor for BaseCounterRow object.
'''
cdef class QualityCounterRow
cdef QualityCounterRow makeQualityCounterRow(
                                          char * ref,
                                          int position,
                                          int num_reads,
                                          char * bases,
                                          double * base_quals,
                                          double * map_quals
                                          ):
     cdef QualityCounterRow row = QualityCounterRow.__new__(QualityCounterRow)
    
     row._ref = ref
     row._position = position
     row._num_reads = num_reads
     row._bases = bases
     row._base_quals = base_quals
     row._map_quals = map_quals     
     
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
