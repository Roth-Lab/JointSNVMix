'''
Created on 2012-01-26

@author: Andrew Roth
'''
#=======================================================================================================================
# Row Factories
#=======================================================================================================================
cdef class RowFactory(object):
    def __init__(self, FastaFile ref_genome, int min_base_qual, int min_map_qual, bint qualities):
        self._min_base_qual = min_base_qual
        self._min_map_qual = min_map_qual
        
        self._ref_genome = ref_genome        
        
        if qualities:
            self._data_factory = QualityDataFactory(min_base_qual, min_map_qual)
        else:
            self._data_factory = CountDataFactory(min_base_qual, min_map_qual)

    cdef JointBinaryCounterRow get_row(self, char * ref, int pos,
                                       PileupColumn normal_column, PileupColumn tumour_column):
                                       
        cdef JointBinaryCounterRow row = JointBinaryCounterRow.__new__(JointBinaryCounterRow)
        
        row._ref = ref
        row._pos = pos
        
        row._ref_base = self._ref_genome.get_reference_base(ref, pos)       
        
        row._var_base = self._get_var_base(row._ref_base, normal_column, tumour_column)
        
        row._data = self._data_factory.get_data(row._ref_base, row._var_base, normal_column, tumour_column)
         
        return row
    
    cdef char * _get_var_base(self, char * ref_base, PileupColumn normal_column, PileupColumn tumour_column):
        '''
        Sort the bases by number observed and return the most common non-reference base.
        
        Return N if the counts of most common non-reference is 0.
        ''' 
        cdef int non_ref_index, normal_counts, tumour_counts   
        cdef base_counts_struct[4] counts
    
        counts[0].base = "A"
        counts[1].base = "C"
        counts[2].base = "G"
        counts[3].base = "T"
        
        for i in range(4):
            normal_counts = normal_column.get_nucleotide_count(counts[i].base,
                                                               self._min_base_qual,
                                                               self._min_map_qual)
            
            tumour_counts = tumour_column.get_nucleotide_count(counts[i].base,
                                                               self._min_base_qual,
                                                               self._min_map_qual) 
            
            counts[i].counts = normal_counts + tumour_counts
    
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

#=======================================================================================================================
# Data factories
#=======================================================================================================================
cdef class DataFactory(object):
    def __init__(self, int min_base_qual, int min_map_qual):
        self._min_base_qual = min_base_qual
        self._min_map_qual = min_map_qual

    cdef JointBinaryData get_data(self, char * ref_base, char * var_base,
                                  PileupColumn normal_column, PileupColumn tumour_column):
        pass       

cdef class CountDataFactory(DataFactory):
    cdef JointBinaryData get_data(self, char * ref_base, char * var_base,
                                  PileupColumn normal_column, PileupColumn tumour_column):
        cdef JointBinaryCountData data = JointBinaryCountData.__new__(JointBinaryCountData)

        data._a_N = normal_column.get_nucleotide_count(ref_base, self._min_base_qual, self._min_map_qual)
        data._a_T = tumour_column.get_nucleotide_count(ref_base, self._min_base_qual, self._min_map_qual)
        
        if strcmp(var_base, 'N') == 0:
            data._b_N = 0
            data._b_T = 0
        else: 
            data._b_N = normal_column.get_nucleotide_count(var_base, self._min_base_qual, self._min_map_qual)
            data._b_T = tumour_column.get_nucleotide_count(var_base, self._min_base_qual, self._min_map_qual)

        return data

cdef class QualityDataFactory(DataFactory):
    cdef JointBinaryData get_data(self, char * ref_base, char * var_base,
                                  PileupColumn normal_column, PileupColumn tumour_column):
        cdef int d_N, d_T
        
        cdef JointBinaryQualityData data = JointBinaryQualityData.__new__(JointBinaryQualityData)

        # Get the number of ref and non-ref bases in tumour.        
        data._a_N = normal_column.get_nucleotide_count(ref_base, self._min_base_qual, self._min_map_qual)
        data._a_T = tumour_column.get_nucleotide_count(ref_base, self._min_base_qual, self._min_map_qual)
        
        if strcmp(var_base, 'N') == 0:                
            data._b_N = 0
            data._b_T = 0
        else:
            data._b_N = normal_column.get_nucleotide_count(var_base, self._min_base_qual, self._min_map_qual)                
            data._b_T = tumour_column.get_nucleotide_count(var_base, self._min_base_qual, self._min_map_qual)
            
        data._d_N = data._a_N + data._b_N
        data._d_T = data._a_T + data._b_T
        
        # Initialise arrays to store base and mapping probabilities. 
        data._q_N = < double *> malloc(sizeof(double) * data._d_N)
        data._r_N = < double *> malloc(sizeof(double) * data._d_N)
        
        data._q_T = < double *> malloc(sizeof(double) * data._d_T)
        data._r_T = < double *> malloc(sizeof(double) * data._d_T)
        
        # Load aligment probabilities.
        self._get_aligment_probabilities(ref_base, var_base, data._q_N, data._r_N, normal_column)
        
        self._get_aligment_probabilities(ref_base, var_base, data._q_T, data._r_T, tumour_column)

        return data
    
    cdef _get_aligment_probabilities(self,
                                     char * ref_base,
                                     char * var_base,
                                     double * q,
                                     double * r,
                                     PileupColumn column):
        '''
        Extract mapping and base qualities for reference position from pileup column.
        
        Since we are only given probabilities for observed base, we assume the remaining probability is evenly shared among
        the 3 other bases.
        
        This function ignores positions which do not match the reference or variant base.
        '''
        cdef char base_char, ref_base_char, var_base_char
        cdef int read_index, i, bq, mq
        cdef double prob
        cdef bint var_base_is_not_N
        
        i = 0
        
        ref_base_char = ref_base[0]
        var_base_char = var_base[0]
        
        if strcmp(var_base, 'N') == 0:
            var_base_is_not_N = 0
        else:
            var_base_is_not_N = 1
        
        for read_index in range(column._depth):
            bq = column._base_quals[read_index]
            mq = column._map_quals[read_index]
                    
            if (bq < self._min_base_qual) or (mq < self._min_map_qual):
                continue
            
            base_char = column._bases[read_index]
    
            if ref_base_char == base_char:
                prob = convert_phred_qual_to_prob(bq)
                q[i] = prob
                
                r[i] = convert_phred_qual_to_prob(mq)
                
                i += 1
            
            elif var_base_char == base_char and var_base_is_not_N:
                prob = convert_phred_qual_to_prob(bq)
                q[i] = (1 - prob) / 3
                
                r[i] = convert_phred_qual_to_prob(mq)
                
                i += 1    

#=======================================================================================================================
# Row objects
#=======================================================================================================================
cdef class JointBinaryCounterRow(object):
    def __dealloc__(self):
        free(self._ref_base)

    def __str__(self):
        out_row = [
                   self.ref,
                   str(self.position),
                   self.ref_base,
                   self.var_base,
                   str(self.normal_ref_counts),
                   str(self.normal_var_counts),
                   str(self.tumour_ref_counts),
                   str(self.tumour_var_counts)
                   ]
        
        
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
        
    property normal_depth:
        def __get__(self):
            return self._data.normal_depth

    property tumour_depth:
        def __get__(self):
            return self._data.tumour_depth
    
    property normal_ref_counts:
        def __get__(self):
            return self._data.normal_ref_counts

    property normal_var_counts:
        def __get__(self):
            return self._data.normal_var_counts        

    property tumour_ref_counts:
        def __get__(self):
            return self._data.tumour_ref_counts

    property tumour_var_counts:
        def __get__(self):
            return self._data.tumour_var_counts
        
    property counts:
        def __get__(self):
            return (self.normal_ref_counts, self.normal_var_counts, self.tumour_ref_counts, self.tumour_var_counts)
    
    property data:
        def __get__(self):
            return self._data 
        
#=======================================================================================================================
# Data object 
#=======================================================================================================================
cdef class JointBinaryData(object):
    property normal_depth:
        def __get__(self):
            return self._a_N + self._b_N

    property tumour_depth:
        def __get__(self):
            return self._a_T + self._b_T
    
    property normal_ref_counts:
        def __get__(self):
            return self._a_N

    property normal_var_counts:
        def __get__(self):
            return self._b_N
        
    property tumour_ref_counts:
        def __get__(self):
            return self._a_T

    property tumour_var_counts:
        def __get__(self):
            return self._b_T

cdef class JointBinaryCountData(JointBinaryData):
    def __init__(self, a_N, b_N, a_T, b_T):
        self._a_N = a_N
        self._b_N = b_N
        self._a_T = a_T
        self._b_T = b_T          

cdef class JointBinaryQualityData(JointBinaryData):
    def __dealloc__(self):
        free(self._q_N)
        free(self._r_N)
        free(self._q_T)
        free(self._r_T)
    
    property normal_base_qualities:
        def __get__(self):
            return [x for x in self._q_N[:self._d_N]]
    
    property tumour_base_qualities:
        def __get__(self):
            return [x for x in self._q_T[:self._d_T]]
        
    property normal_mapping_qualities:
        def __get__(self):
            return [x for x in self._r_N[:self._d_N]]
        
    property tumour_mapping_qualities:
        def __get__(self):
            return [x for x in self._r_T[:self._d_T]]
        
#===============================================================================
# Utility functions for finding non-ref bases and counts.
#===============================================================================
cdef int compare_base_counts_struct(const_void * a, const_void * b):
    cdef base_counts_struct * first
    cdef base_counts_struct * second
    
    first = < base_counts_struct *> a
    second = < base_counts_struct *> b

    return first.counts - second.counts

cpdef double convert_phred_qual_to_prob(int qual):
    '''
    Converts integer mapping or base quality on the phred scale to a probability.
    '''
    cdef double base, exponent, prob
    
    exponent = -1 * (< double > qual) / 10
    base = 10
    
    prob = 1 - pow(base, exponent)

    return prob