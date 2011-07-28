cdef class ThresholdClassifier(Classifier):    
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return ThresholdClassifierRefIterator(
                                              ref,
                                              self._counter.iter_ref(ref),
                                              **kwargs
                                              )
             
cdef class ThresholdClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryBaseCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        self._normal_threshold = kwargs.get('normal_threshold', 0.05)
        self._tumour_threshold = kwargs.get('tumour_threshold', 0.1)

    cdef tuple _get_labels(self):
        cdef int normal_genotype, tumour_genotype, joint_genotype 
        cdef list labels
        cdef JointBinaryCounterRow row
        
        row = self._iter._current_row
        
        normal_genotype = self._get_genotype(row._normal_counts.A, row._normal_counts.B, self._normal_threshold)
        tumour_genotype = self._get_genotype(row._tumour_counts.A, row._tumour_counts.B, self._tumour_threshold)        
        
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
