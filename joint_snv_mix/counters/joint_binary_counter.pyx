'''
Use the convention every position at C-level is 0-based. If the position is
accessed by property mechanism from python layer it should be one based.

Created on 2011-06-29

@author: Andrew Roth
'''
cdef class JointBinaryBaseCounter(Counter):
    '''
    Class for iterating over positions from paired genome files counting bases.
    '''
    def __init__(self, Samfile normal_bam, Samfile tumour_bam, Fastafile ref_genome_fasta, min_base_qual=10, min_map_qual=10):
        self._normal_counter = BaseCounter(normal_bam, min_base_qual, min_map_qual)
        self._tumour_counter = BaseCounter(tumour_bam, min_base_qual, min_map_qual)
        
        self._ref_genome_fasta = ref_genome_fasta
        
        self._refs = tuple(set(self._normal_counter.refs) & set(self._tumour_counter.refs)) 
        
    def iter_ref(self, ref):
        if ref not in self.refs:
            raise Exception("Invalid reference passed.")
        
        return JointBinaryBaseCounterIterator(
                                               ref,
                                               self._normal_counter.iter_ref(ref),
                                               self._tumour_counter.iter_ref(ref),
                                               self._ref_genome_fasta
                                               )
        
cdef class JointBinaryBaseCounterIterator(JointRefIterator):
    def __init__(self,
                 char * ref,
                 BaseCounterRefIterator normal_iter,
                 BaseCounterRefIterator tumour_iter,
                 Fastafile ref_genome_fasta):
        
        self._ref = ref
        
        self._normal_iter = normal_iter
        self._tumour_iter = tumour_iter
        
        self._ref_genome_fasta = ref_genome_fasta
        
        self._position = -1
    
    cdef parse_current_position(self):
        cdef int region_length
        cdef char * ref_base
        cdef BaseCounterRow normal_row
        cdef BaseCounterRow tumour_row
        
        self._normal_iter.parse_current_position()        
        normal_row = self._normal_iter._current_row
        
        self._tumour_iter.parse_current_position()
        tumour_row = self._tumour_iter._current_row
        
        region_length = 1 
        
        ref_base = self._ref_genome_fasta._fetch(
                                                 normal_row._ref,
                                                 normal_row._position,
                                                 normal_row._position + 1,
                                                 & region_length
                                                 )
    
        self._current_row = makeJointBinaryCounterRow(ref_base, normal_row, tumour_row)
    
cdef class JointBinaryCounterRow
cdef JointBinaryCounterRow makeJointBinaryCounterRow(char * ref_base, BaseCounterRow normal_row, BaseCounterRow tumour_row):
    '''
    Constructor method for creating a JointBinaryCounterRow from C.
    '''
    cdef char * non_ref_base

    cdef JointBinaryCounterRow row = JointBinaryCounterRow.__new__(JointBinaryCounterRow)
    
    row._ref = normal_row._ref
    row._position = normal_row._position
    
    row._ref_base = ref_base
    
    non_ref_base = get_non_ref_base(ref_base, normal_row, tumour_row)
    
    row._non_ref_base = non_ref_base
     
    row._normal_counts = get_binary_counts(ref_base, non_ref_base, normal_row)
    row._tumour_counts = get_binary_counts(ref_base, non_ref_base, tumour_row)
    
    row._normal_depth = row._normal_counts.A + row._normal_counts.B
    row._tumour_depth = row._tumour_counts.A + row._tumour_counts.B   
     
    return row

cdef class JointBinaryCounterRow(PairedSampleCounterRow):
    '''
    Class for storing binary count data from a pair of Bam files at a position.
    '''    
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")
    
    def __dealloc__(self):
        free(self._ref_base)
    
    def __str__(self):
        '''
        Overide parent class method to include ref and non_ref base.
        '''
        out_row = [self.ref, str(self.position), self.ref_base, self.non_ref_base]
        out_row.extend([str(x) for x in self.counts])
        
        return "\t".join(out_row)
    
    property counts:
        def __get__(self):
            return (
                    self._normal_counts.A,
                    self._normal_counts.B,
                    self._tumour_counts.A,
                    self._tumour_counts.B
                    )
    
    property ref_base:
        def __get__(self):
            return self._ref_base

    property non_ref_base:
        def __get__(self):
            return self._non_ref_base
    
#===============================================================================
# Utility functions for finding non-ref bases and counts.
#===============================================================================
cdef int compare_base_counts_struct(const_void * a, const_void * b):
    cdef base_counts_struct * first
    cdef base_counts_struct * second
    
    first = < base_counts_struct *> a
    second = < base_counts_struct *> b

    return first.counts - second.counts

cdef char * get_non_ref_base(char * ref_base, BaseCounterRow normal_row, BaseCounterRow tumour_row):
    '''
    Sort the bases by number observed and return the most common non-reference base.
    
    Return N if the counts of most common non-reference is 0.
    ''' 
    cdef int non_ref_index    
    cdef base_counts_struct[4] counts

    counts[0].base = "A"
    counts[1].base = "C"
    counts[2].base = "G"
    counts[3].base = "T"
    
    counts[0].counts = normal_row._counts.A + tumour_row._counts.A
    counts[1].counts = normal_row._counts.C + tumour_row._counts.C
    counts[2].counts = normal_row._counts.G + tumour_row._counts.G
    counts[3].counts = normal_row._counts.T + tumour_row._counts.T   

    # Sort the structs by counts field. 
    qsort(counts, 4, sizeof(base_counts_struct), compare_base_counts_struct)
            
    # If the most prevalent base is not the reference return it.
    if strcmp(counts[3].base, ref_base) != 0:
        non_ref_index = 3
    else:
        non_ref_index = 2
    
    if counts[non_ref_index].counts > 0:
        return counts[non_ref_index].base
    else:
        return 'N'
        
cdef binary_counts_struct get_binary_counts(char * ref_base, char * non_ref_base, BaseCounterRow row):
    cdef binary_counts_struct counts

    counts.A = get_counts(ref_base, row._counts)
    counts.B = get_counts(non_ref_base, row._counts)
    
    return counts

cdef int get_counts(char * base, counts_struct counts):
    cdef int result

    if strcmp(base, "A") == 0:
        result = counts.A        
    elif strcmp(base, "C") == 0:
        result = counts.C
    elif strcmp(base, "G") == 0:
        result = counts.G
    elif strcmp(base, "T") == 0:
        result = counts.T
    
    return result
