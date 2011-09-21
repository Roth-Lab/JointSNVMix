cdef class SnvMixData(object):
    pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneData(SnvMixData):
    pass

cdef SnvMixOneData makeSnvMixOneData(binary_counts_struct counts):
    cdef SnvMixOneData data = SnvMixOneData.__new__(SnvMixOneData)
    
    data.counts[0] = counts.A
    data.counts[1] = counts.B
    
    return data

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoData(SnvMixData):
    def __dealloc__(self):
        free(self.q)
        free(self.r)

cdef SnvMixTwoData makeSnvMixTwoData(base_map_qualities_struct data_struct):
    cdef read_index, i, l
    cdef double temp_q
    cdef SnvMixTwoData data = SnvMixTwoData.__new__(SnvMixTwoData)
    
    l = data_struct.depth.A + data_struct.depth.B
    
    data.depth = l
    
    data.q = < double *> malloc(l * sizeof(double))
    data.r = < double *> malloc(l * sizeof(double))
    
    i = 0
    
    for read_index in range(data_struct.depth.A):
        data.q[i] = get_phred_qual_to_prob(data_struct.base_quals.A[read_index])
        data.r[i] = get_phred_qual_to_prob(data_struct.map_quals.A[read_index])
        
        i += 1
    
    for read_index in range(data_struct.depth.B):        
        temp_q = get_phred_qual_to_prob(data_struct.base_quals.B[read_index])
        data.q[i] = (1 - temp_q) / 3
                
        data.r[i] = get_phred_qual_to_prob(data_struct.map_quals.B[read_index])
        
        i += 1
    
    return data

cdef double get_phred_qual_to_prob(int qual):
    cdef double base, exp, prob
    
    exp = -1 * (< double > qual) / 10
    base = 10
    
    prob = 1 - pow(base, exp)

    return prob
