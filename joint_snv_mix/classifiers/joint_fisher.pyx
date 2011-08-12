cdef class JointFisherClassifier(Classifier):    
    def __init__(self, ** kwargs):
        self._p_value_threshold = kwargs.get('p_value_threshold', 0.05)
        
        self._indep_classifier = IndependentFisherClassifier(**kwargs)

    cdef tuple _get_labels(self, PairedSampleBinomialCounterRow row):
        cdef PValues pv
        cdef tuple indep_labels
        cdef list corrected_labels
        
        indep_labels = self._indep_classifier._get_labels(row)
        
        # Correct somatic calls by using fisher exact test to see if normal/tumour counts are significantly different.
        if indep_labels[1] == 1 or indep_labels[2] == 1:
            pv = fisher_exact_test(
                                   (< JointBinaryCounterRow > row)._normal_counts.A,
                                   (< JointBinaryCounterRow > row)._normal_counts.B,
                                   (< JointBinaryCounterRow > row)._tumour_counts.A,
                                   (< JointBinaryCounterRow > row)._tumour_counts.B
                                   )

            # Counts are not significantly different
            if pv.two_tail > self._p_value_threshold:            
                corrected_labels = [0] * 9
                corrected_labels[0] = 1
                
                indep_labels = tuple(corrected_labels)

        return indep_labels
