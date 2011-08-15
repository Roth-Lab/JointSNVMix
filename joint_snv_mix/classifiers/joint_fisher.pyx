DEF NUM_JOINT_GENOTYPES = 9

cdef class JointFisherClassifier(Classifier): 
    def __init__(self, ** kwargs):
        self._p_value_threshold = kwargs.get('p_value_threshold', 0.05)
        
        self._indep_classifier = IndependentFisherClassifier(**kwargs)

    cdef double * _get_labels(self, PairedSampleBinomialCounterRow row):
        cdef int g
        cdef PValues pv
        cdef double * labels
        
        labels = self._indep_classifier._get_labels(row)
        
        # Correct somatic calls by using fisher exact test to see if normal/tumour counts are significantly different.
        if labels[1] == 1 or labels[2] == 1:
            pv = fisher_exact_test(
                                   row._counts[0],
                                   row._counts[1],
                                   row._counts[2],
                                   row._counts[3]
                                   )

            # Counts are not significantly different
            if pv.two_tail > self._p_value_threshold:
                for g in range(NUM_JOINT_GENOTYPES):
                    labels[g] = 0
                    
                labels[0] = 1
                
        return labels
