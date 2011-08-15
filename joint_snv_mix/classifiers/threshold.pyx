DEF NUM_JOINT_GENOTYPES = 9

cdef class ThresholdClassifier(Classifier):
    def __init__(self, **kwargs):
        self._normal_threshold = kwargs.get('normal_threshold', 0.05)
        self._tumour_threshold = kwargs.get('tumour_threshold', 0.05)

    cdef double * _get_labels(self, PairedSampleBinomialCounterRow row):
        cdef int normal_genotype, tumour_genotype, joint_genotype, g 
        cdef double * labels
        
        normal_genotype = self._get_genotype(row._counts[0],
                                             row._counts[1],
                                             self._normal_threshold)
        
        tumour_genotype = self._get_genotype(row._counts[2],
                                             row._counts[3],
                                             self._tumour_threshold)        
        
        joint_genotype = 3 * normal_genotype + tumour_genotype
        
        labels = < double *> malloc(NUM_JOINT_GENOTYPES * sizeof(double))
        
        for g in range(NUM_JOINT_GENOTYPES):
            if g == joint_genotype:
                labels[g] = 1
            else:
                labels[g] = 0
        
        return labels
                
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
