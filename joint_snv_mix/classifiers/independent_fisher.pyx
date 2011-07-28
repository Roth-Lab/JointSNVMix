cdef class IndependentFisherClassifier(Classifier):    
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return IndependentFisherClassifierRefIterator(
                                                      ref,
                                                      self._counter.iter_ref(ref),
                                                      **kwargs
                                                      )
             
cdef class IndependentFisherClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryBaseCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        self._min_var_freq = kwargs.get('min_var_freq', 0.1)
        self._hom_var_freq = kwargs.get('hom_var_freq', 0.9)
        self._p_value_threshold = kwargs.get('p_value_threshold', 0.05)
        self._expected_error_rate = kwargs.get('expected_error_rate', 0.001)

    cdef tuple _get_labels(self):
        cdef int normal_genotype, tumour_genotype, joint_genotype 
        cdef list labels
        cdef JointBinaryCounterRow row
        
        row = self._iter._current_row
        
        normal_genotype = self._get_genotype(row._normal_counts.A, row._normal_counts.B)
        tumour_genotype = self._get_genotype(row._tumour_counts.A, row._tumour_counts.B)        
        
        joint_genotype = 3 * normal_genotype + tumour_genotype
        
        labels = [0] * 9
        
        labels[joint_genotype] = 1        
        
        return tuple(labels)
                
    cdef int _get_genotype(self, int a, int b):
        '''
        Given count data call genotype. 0 - No variant, 1 - Het variant, 2 - Homozygous variant.
        '''
        cdef int d
        cdef float pv, var_freq
        
        d = a + b
        
        pv = self._get_significance(a, b)
        
        if d > 0:
            var_freq = (< float > b) / d
        else:
            var_freq = 0
        
        # Counts differ from expectation due to random error and depth of coverage is acceptable.
        if pv <= self._p_value_threshold:
            # Homozygous variant.
            if var_freq >= self._hom_var_freq:
                return 2
            # Heterozygous variant.
            elif var_freq >= self._min_var_freq:
                return 1
        
        # Something failed so call non-variant.        
        
        return 0
        
    cdef float _get_significance(self, int a, int b):
        cdef int expected_a, expected_b, d
        cdef PValues pv
        
        d = a + b
        
        expected_b = < int > floor(self._expected_error_rate * d)
        expected_a = d - expected_b
        
        pv = fisher_exact_test(
                               expected_a,
                               expected_b,
                               a,
                               b
                               )        
        return pv.right_tail
