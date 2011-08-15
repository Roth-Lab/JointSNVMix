'''
Use the convention every position at C-level is 0-based. If the position is
accessed by property mechanism from python layer it should be one based.

Created on 2011-07-28

@author: Andrew Roth
'''
cdef class JointBinaryQualityCounter(Counter):
    '''
    Class for iterating over positions from paired genome files counting bases.
    '''
    def __init__(self, Samfile normal_bam, Samfile tumour_bam, Fastafile ref_genome_fasta):
        self._normal_counter = QualityCounter(normal_bam)
        self._tumour_counter = QualityCounter(tumour_bam)
        
        self._ref_genome_fasta = ref_genome_fasta
        
        self._refs = tuple(set(self._normal_counter.refs) & set(self._tumour_counter.refs)) 
        
    def iter_ref(self, ref):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return JointBinaryQualityCounterIterator(
                                                 ref,
                                                 self._normal_counter.iter_ref(ref),
                                                 self._tumour_counter.iter_ref(ref),
                                                 self._ref_genome_fasta
                                                 )
        
cdef class JointBinaryQualityCounterIterator(JointRefIterator):
    def __init__(self,
                 char * ref,
                 QualityCounterRefIterator normal_iter,
                 QualityCounterRefIterator tumour_iter,
                 Fastafile ref_genome_fasta):
        
        self._ref = ref
        
        self._normal_iter = normal_iter
        self._tumour_iter = tumour_iter
        
        self._ref_genome_fasta = ref_genome_fasta
        
        self._position = -1

    cdef parse_current_position(self):
        cdef int region_length
        cdef char * ref_base
        cdef QualityCounterRow normal_row
        cdef QualityCounterRow tumour_row
        
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
    
        self._current_row = makeJointBinaryQualityCounterRow(ref_base, normal_row, tumour_row)

cdef class JointBinaryQualityCounterRow(PairedSampleBinomialCounterRow):
    '''
    Class for storing binary count data with qualities from a pair of Bam files at a position.
    '''    
    def __dealloc__(self):
        destroy_base_map_qualities_struct(self._normal_data)
        destroy_base_map_qualities_struct(self._tumour_data)
        free(self._ref_base)

    property counts:
        def __get__(self):
            return (
                    self._normal_data.depth.A,
                    self._normal_data.depth.B,
                    self._tumour_data.depth.A,
                    self._tumour_data.depth.B,
                    )

#=======================================================================================================================
# Row factory function and helper functions.
#=======================================================================================================================
cdef JointBinaryQualityCounterRow makeJointBinaryQualityCounterRow(char * ref_base,
                                                                   QualityCounterRow normal_row,
                                                                   QualityCounterRow tumour_row):
    '''
    Constructor method for creating a JointBinaryQualityCounterRow from C.
    '''
    cdef char * non_ref_base

    cdef JointBinaryQualityCounterRow row = JointBinaryQualityCounterRow.__new__(JointBinaryQualityCounterRow)
    
    row._ref = normal_row._ref
    row._position = normal_row._position
    
    row._ref_base = ref_base
    
    non_ref_base = get_non_ref_base(ref_base, normal_row, tumour_row)
    
    row._non_ref_base = non_ref_base
     
    row._normal_data = get_base_map_qualities(ref_base, non_ref_base, normal_row)
    row._tumour_data = get_base_map_qualities(ref_base, non_ref_base, tumour_row)
    
    row._normal_depth = row._normal_data.depth.A + row._normal_data.depth.B
    row._tumour_depth = row._tumour_data.depth.A + row._tumour_data.depth.B
    
    row._counts[0] = row._normal_data.depth.A
    row._counts[1] = row._normal_data.depth.B
    row._counts[2] = row._tumour_data.depth.A
    row._counts[3] = row._tumour_data.depth.B
     
    return row

cdef base_map_qualities_struct get_base_map_qualities(char * ref_base, char * non_ref_base, QualityCounterRow row):
    cdef int i, ref_counts, non_ref_counts, A_index, B_index
    cdef base_map_qualities_struct data
    
    ref_counts = row.get_counts(ref_base)
    non_ref_counts = row.get_counts(non_ref_base)
    
    data = create_base_map_qualities_struct(ref_counts, non_ref_counts)
    
    A_index = 0 
    B_index = 0
    
    for i in range(row._depth):
        if row._bases[i] == ref_base[0]:
            data.base_quals.A[A_index] = row._base_quals[i]
            
            data.map_quals.A[A_index] = row._map_quals[i]
            
            A_index += 1
            
        elif row._bases[i] == non_ref_base[0]:
            data.base_quals.B[B_index] = row._base_quals[i]
            
            data.map_quals.B[B_index] = row._map_quals[i]
            
            B_index += 1
    
    return data

cdef base_map_qualities_struct create_base_map_qualities_struct(int A_counts, int B_counts):
    cdef base_map_qualities_struct data
    
    data.depth.A = A_counts
    data.depth.B = B_counts

    data.base_quals.A = < int * > malloc(A_counts * sizeof(int)) 
    data.base_quals.B = < int * > malloc(B_counts * sizeof(int))
    
    data.map_quals.A = < int * > malloc(A_counts * sizeof(int))
    data.map_quals.B = < int * > malloc(B_counts * sizeof(int))
    
    return data

cdef void destroy_base_map_qualities_struct(base_map_qualities_struct data):
    free(data.base_quals.A)
    free(data.base_quals.B)
    
    free(data.map_quals.A)
    free(data.map_quals.B)
