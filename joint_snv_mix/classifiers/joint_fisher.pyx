cdef class JointFisherClassifier(Classifier):    
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return JointFisherClassifierRefIterator(
                                                ref,
                                                self._counter.iter_ref(ref),
                                                **kwargs
                                                )
             
cdef class JointFisherClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryBaseCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        self._min_var_freq = kwargs.get('min_var_freq', 0.1)
        self._hom_var_freq = kwargs.get('hom_var_freq', 0.9)
        self._p_value_threshold = kwargs.get('p_value_threshold', 0.05)
        self._expected_error_rate = kwargs.get('expected_error_rate', 0.001)
        
        self._indep_iter = IndependentFisherClassifierRefIterator(
                                                                  ref,
                                                                  iter,
                                                                  **kwargs
                                                                  )

    cdef tuple _get_labels(self):
        cdef PValues pv
        cdef tuple indep_labels
        cdef list corrected_labels
        cdef JointBinaryCounterRow row
        
        indep_labels = self._indep_iter._get_labels()
        
        row = self._iter._current_row
        
        # Correct somatic calls by using fisher exact test to see if normal/tumour counts are significantly different.
        if indep_labels[1] == 1 or indep_labels[2] == 1:
            pv = fisher_exact_test(
                                   row._normal_counts.A,
                                   row._normal_counts.B,
                                   row._tumour_counts.A,
                                   row._tumour_counts.B
                                   )

            # Counts are not significantly different
            if pv.two_tail > self._p_value_threshold:            
                corrected_labels = [0] * 9
                corrected_labels[0] = 1
                
                return tuple(corrected_labels)
                
                
        
        return indep_labels
