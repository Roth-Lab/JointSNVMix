cdef class ThresholdClassifier(Classifier):    
    def __init__(self, **kwargs):
        self._normal_threshold = kwargs.get('normal_threshold', 0.05)
        self._tumour_threshold = kwargs.get('tumour_threshold', 0.05)

    cdef tuple _get_labels(self):
        cdef int normal_genotype, tumour_genotype, joint_genotype 
        cdef list labels
        
        normal_genotype = self._get_genotype((< JointBinaryCounterRow > row)._normal_counts.A, 
                                             (< JointBinaryCounterRow > row)._normal_counts.B, 
                                             self._normal_threshold)
        
        tumour_genotype = self._get_genotype((< JointBinaryCounterRow > row)._tumour_counts.A, 
                                             (< JointBinaryCounterRow > row)._tumour_counts.B, 
                                             self._tumour_threshold)        
        
        joint_genotype = 3 * normal_genotype + tumour_genotype
        
        labels = [0] * 9
        
        labels[joint_genotype] = 1        
        
        return tuple(labels)
                
    cdef int _get_genotype(self, int a, int b, float freq_threshold):
        '''
        Given count data call genotype. 0 - No variant, 1 - Het variant, 2 - Homozygous variant.
        '''
        cdef int d
        cdef float var_freq
        
        d = a + b
        
        if d > 0:
            var_freq = (< float > b) / d
        else:
            var_freq = 0
        
        # Homozygous variant.
        if var_freq >= freq_threshold:
            if var_freq >= (1 - freq_threshold):
                return 2
            else:
                return 1
        else:
            return 0
