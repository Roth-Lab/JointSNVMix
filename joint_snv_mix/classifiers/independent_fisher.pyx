cdef class IndependentFisherClassifier(Classifier):    
    def __init__(self, **kwargs):        
        self._min_var_freq = kwargs.get('min_var_freq', 0.1)
        self._hom_var_freq = kwargs.get('hom_var_freq', 0.9)
        self._p_value_threshold = kwargs.get('p_value_threshold', 0.05)
        self._expected_error_rate = kwargs.get('expected_error_rate', 0.001)
        self._min_var_depth = kwargs.get('min_var_depth', 4)

    cdef tuple _get_labels(self, PairedSampleBinomialCounterRow row):
        cdef int normal_genotype, tumour_genotype, joint_genotype 
        cdef list labels
        
        normal_genotype = self._get_genotype((< JointBinaryCounterRow > row)._normal_counts.A,
                                             (< JointBinaryCounterRow > row)._normal_counts.B)
        
        tumour_genotype = self._get_genotype((< JointBinaryCounterRow > row)._tumour_counts.A,
                                             (< JointBinaryCounterRow > row)._tumour_counts.B)        
        
        joint_genotype = 3 * normal_genotype + tumour_genotype
        
        labels = [0] * 9
        
        labels[joint_genotype] = 1        
        
        return tuple(labels)
                
    cdef int _get_genotype(self, int a, int b):
        '''
        Given count data call genotype. 0 - No variant, 1 - Het variant, 2 - Homozygous variant.
        '''
        cdef int d, genotype
        cdef float pv, var_freq,
        
        d = a + b
        
        pv = self._get_significance(a, b)
        
        if d > 0:
            var_freq = (< float > b) / d
        else:
            var_freq = 0
        
        # By default assume AA genotype.
        genotype = 0
        
        # Check if counts differ from expectation due to random error and depth of coverage is acceptable.
        if pv <= self._p_value_threshold:
            # Homozygous BB variant.
            if var_freq >= self._hom_var_freq:
                genotype = 2
            # Heterozygous AB variant.
            elif var_freq >= self._min_var_freq:
                genotype = 1
        
        # Final check to make sure there are enough variant reads otherwise call AA.
        if b < self._min_var_depth:
            genotype = 0
        
        # Something failed so call non-variant.        
        return genotype
        
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
