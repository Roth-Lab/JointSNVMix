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
    
cdef class JointBinaryCounterRow(PairedSampleCounterRow):
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
            
#=======================================================================================================================
# Row factory function
#=======================================================================================================================
cdef JointBinaryCounterRow makeJointBinaryCounterRow(char * ref_base,
                                                     BaseCounterRow normal_row,
                                                     BaseCounterRow tumour_row):
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
