'''
Created on 2012-01-26

@author: Andrew Roth
'''
cdef class FeatureExtractor(object):
    def __init__(self, BamFile normal_bam, BamFile tumour_bam):
        self._normal_bam = normal_bam
        self._tumour_bam = tumour_bam
    
    def get_features(self, JointBinaryCounterRow row):
        normal_features = self._get_genome_features(row, self._normal_bam)
        tumour_features = self._get_genome_features(row, self._tumour_bam)
        
        return normal_features + tumour_features
    
    cdef list _get_genome_features(self, JointBinaryCounterRow row, BamFile bam_file):
        pileup_iter = bam_file.get_pileup_iterator(row._ref, start=row._pos, stop=row._pos + 1)
        
        while True:
            pileup_iter.cnext()
            
            if pileup_iter._pos = row._pos:
                break
        
        pileup_column = pileup_iter.get_extended_pileup_column()
        
        features = self._parse_column(row, pileup_column)
        
        return features
    
    cdef list _parse_column(self, JointBinaryCounterRow row, ExtendedPileupColumn pileup_column):
        depth = pileup_column._depth
        
        ref_base = row._ref_base
        var_base = row._var_base
        
        ref_base_features = self._get_base_features(ref_base, pileup_column)
        var_base_features = self._get_base_features(var_base, pileup_column)
        
        snv_mix_probs = self._model.predict(row)
        
        return ref_base_feature + var_base_features + snv_mix_probs
    
    cdef list _get_base_features(self, char * base, ExtendedPileupColumn pileup_column):
        base_quals_sum = self._get_base_quals_sum(base, pileup_column, 1)        
        base_quals_square_sum = self._get_base_quals_sum(base, pileup_column, 2)
        
        map_quals_sum = self._get_map_quals_sum(ref_base, pileup_column, 1)
        map_quals_square_sum = self._get_map_quals_sum(ref_base, pileup_column, 2)
        
        forward_strand_count = self._get_forward_strand_count(base, pileup_column)                
        reverse_strand_count = pileup_column.get_nucleotide_count(base, 0, 0) - forward_strand_count 
        
        tail_distance_sum = self._get_tail_distance_sum(base, pileup_column, 1)
        tail_distance_square_sum = self._get_tail_distance_sum(base, pileup_column, 2)
                
        return [
                base_quals_sum,
                base_quals_square_sum,
                map_quals_sum,
                map_quals_square_sum,
                forward_strand_count,
                reverse_strand_count,
                tail_distance_sum,
                tail_distance_square_sum
                ]
        
    cdef double _get_base_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double base_qual_sum
        
        base_qual_sum = 0
    
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._base[i]:
                base_qual_sum += pow(pileup_column._base_quals[i], exponent)
        
        return base_qual_sum

    cdef double _get_map_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double map_qual_sum
        
        map_qual_sum = 0
    
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._base[i]:
                map_qual_sum += pow(pileup_column._map_quals[i], exponent)
        
        return map_qual_sum
    
    cdef int _get_forward_strand_count(self, char * base, ExtendedPileupColumn pileup_column):
        cdef int i, forward_strand_count
        
        forward_strand_count = 0
        
        for i in range(pileup_column._depth):
            if (base[0] == pileup_column._base[i]) and pileup_column._is_forward_strand[i]:
                forward_strand_count += 1
        
        return forward_strand_count
    
    cdef double _get_tail_distance_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double tail_distance_sum
        
        tail_distance_sum = 0
        
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._base[i]:
                tail_distance_sum += pow(pileup_column._tail_distance[i], exponent)
        
        return tail_distance_sum
                 
    
            
        
    
