'''
Classes for creating iterators for count data over a pair of genomes. 

Created on 2012-01-18

@author: Andrew Roth
'''
cdef class JointBinaryCounter(Counter):
    '''
    Class for iterating over positions from paired genome files counting bases.
    
    Parameters :
        type - What type of counter. Options : base, quality
        normal_counter - A Counter object for the normal genome.
        tumour_counter - A Counter object for the tumour genome.
        ref_genome - A FastaFile for the reference genome.
    '''
    def __init__(self,
                 char * type,
                 BamFile normal_bam,
                 BamFile tumour_bam,
                 FastaFile ref_genome,
                 int min_base_qual=0,
                 int min_map_qual=0):
        self._type = type
        
        self._normal_bam = normal_bam
        self._tumour_bam = tumour_bam
        
        self._ref_genome = ref_genome
        
        self._refs = tuple(set(self._normal_bam.refs) & set(self._tumour_bam.refs)) 

    property refs:
        '''
        Read only access to list of available references.
        '''
        def __get__(self):
            return self._refs
    
    def get_ref_iterator(self, ref):
        if ref not in self.refs:
            raise Exception("Invalid reference passed.")
        
        return JointBinaryCounterIterator(self._type,
                                          ref,
                                          self._normal_bam.get_pileup_iterator(ref),
                                          self._tumour_bam.get_pileup_iterator(ref),
                                          self._ref_genome,
                                          self._min_base_qual,
                                          self._min_map_qual)
        
cdef class JointBinaryCounterIterator(RefIterator):
    def __init__(self,
                 char * type,
                 char * ref,
                 PileupIterator normal_iter,
                 PileupIterator tumour_iter,
                 FastaFile ref_genome
                 int min_base_qual,
                 int min_map_qual):        
        self._type = type        
        self._ref = ref
        
        self._min_base_qual = 0
        self._min_map_qual = 0
        
        self._normal_iter = normal_iter
        self._tumour_iter = tumour_iter
        
        self._ref_genome = ref_genome
        
        self._pos = -1
        
    def __iter__(self):
        return self
    
    def __next__(self):
        '''
        Python level next() method.
        '''
        self.cnext()
        
        return self._current_row

    property ref:
        '''
        Read only access to reference which the iterator runs over.
        '''
        def __get__(self):
            return self._ref
    
    property position:
        '''
        Read only access to 1-based current position of iterator.
        '''
        def __get__(self):
            return self._position + 1    
        
    cdef advance_position(self):        
        cdef int normal_pos
        cdef int tumour_pos
        
        self._normal_iter.advance_position()        
        self._tumour_iter.advance_position()
                
        while True:
            normal_pos = self._normal_iter._pos
            tumour_pos = self._tumour_iter._pos
            
            if normal_pos == tumour_pos:
                self._pos = normal_pos                                       
                break
            elif normal_pos < tumour_pos:
                self._normal_iter.advance_position()
            elif normal_pos > tumour_pos:
                self._tumour_iter.advance_position()
            else:
                raise Exception("Error in joint pileup iterator.")        
    
    cdef parse_current_position(self):
        cdef char * ref_base
        cdef PileupColumn normal_column
        cdef PileupColumn tumour_column
        
        self._normal_iter.parse_current_position()        
        normal_column = self._normal_iter._current_row
        
        self._tumour_iter.parse_current_position()
        tumour_column = self._tumour_iter._current_row
                
        self._current_row = self._make_counter_row(normal_column, tumour_column)
    
    cdef _make_counter_row(self, normal_column, tumour_column):
        cdef JointBinaryBaseCounterRow row = JointBinaryBaseCounterRow.__new__(JointBinaryBaseCounterRow)
        
        row._ref = self._ref
        row._pos = self._pos
        
        row._ref_base = self._ref_genome.get_reference_base(self._ref, self._pos)       
        row._var_base = get_var_base(ref_base, normal_column, tumour_column, self._min_base_qual, self._min_map_qual)
        
        if self._type == "base":
            row._data = self._make_count_data(row._ref_base, row._var_base, normal_column, tumour_column)
        if self._type == "quality":
            row._data = self._make_quality_data(row._ref_base, normal_column, tumour_column)            
         
        return row

    cdef _make_count_data(self, ref_base, var_base, normal_column, tumour_column):
        cdef JointBinaryCountData data = JointBinaryCountData.__new__(JointBinaryCountData)

        data._a_N = normal_column.get_nucleotide_count(ref_base, self._min_base_qual, self._min_map_qual)
        data._b_N = normal_column.get_nucleotide_count(var_base, self._min_base_qual, self._min_map_qual)

        data._a_T = tumour_column.get_nucleotide_count(ref_base, self._min_base_qual, self._min_map_qual)
        data._b_T = tumour_column.get_nucleotide_count(var_base, self._min_base_qual, self._min_map_qual)
        
        return data
    
    cdef _make_quality_data(self, ref_base, normal_column, tumour_column):
        cdef JointBinaryQualityData data = JointBinaryQualityData.__new__(JointBinaryQualityData)

        data._normal_depth = normal_column._depth
        data._tumour_depth = tumour_column._depth

        data._q_N = < double *> malloc(sizeof(double) * normal_column._depth)
        data._r_N = < double *> malloc(sizeof(double) * normal_column._depth)
        
        data._q_T = < double *> malloc(sizeof(double) * tumour_column._depth)
        data._r_T = < double *> malloc(sizeof(double) * tumour_column._depth) 

        normal_column.get_base_probabilities(ref_base, data._q_N)
        normal_column.get_mapping_probabilities(data._r_N)
        
        tumour_column.get_base_probabilities(ref_base, data._q_T)
        tumour_column.get_mapping_probabilities(data._r_T)        

        return data      
    
#=======================================================================================================================
# Row objects
#=======================================================================================================================
cdef class JointBinaryCounterRow(CounterRow):
    def __dealloc__(self):
        free(self._ref_base)
        free(self._var_base)

    def __str__(self):
        out_row = [self.ref, str(self.position), self.ref_base, self.non_ref_base]
        out_row.extend([str(x) for x in self.counts])
        
        return "\t".join(out_row)    
    
    property ref:
        def __get__(self):
            return self._ref
    
    property position:
        '''
        1-based position
        '''
        def __get__(self):
            return self._pos + 1
    
    property ref_base:
        def __get__(self):
            return self._ref_base
        
    property var_base:
        def __get__(self):
            return self._var_base
    
    property data:
        def __get__(self):
            return self._data 

cdef class JointBinaryBaseCounterRow(JointBinaryCounterRow):
    pass

cdef class JointBinaryQualityCounterRow(JointBinaryCounterRow):
    pass

#=======================================================================================================================
# Data object 
#=======================================================================================================================
cdef class JointBinaryData(object):
    pass

cdef class JointBinaryCountData(object):
    property normal_depth:
        def __get__(self):
            return self._a_N + self._b_N

    property tumour_depth:
        def __get__(self):
            return self._a_T + self._b_T        

cdef class JointBinaryQualityData(object):
    def __dealloc__(self):
        free(self._q_N)
        free(self._r_N)
        free(self._q_T)
        free(self._r_T)
    
    property normal_depth:
        def __get__(self):
            return self._normal_depth
    
    property tumour_depth:
        def __get__(self):
            return self._tumour_depth

#===============================================================================
# Utility functions for finding non-ref bases and counts.
#===============================================================================
cdef int compare_base_counts_struct(const_void * a, const_void * b):
    cdef base_counts_struct * first
    cdef base_counts_struct * second
    
    first = < base_counts_struct *> a
    second = < base_counts_struct *> b

    return first.counts - second.counts

cdef char * get_var_base(char * ref_base,
                         PileupColumn normal_column,
                         PileupColumn tumour_column,
                         int min_base_qual,
                         int min_map_qual):
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
    
    for i in range(4):
        counts[i].counts = normal_column.get_nucleotide_count(counts[i].base, min_base_qual, min_map_qual) + \
                           tumour_column.get_nucleotide_count(counts[i].base, min_base_qual, min_map_qual) 

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

