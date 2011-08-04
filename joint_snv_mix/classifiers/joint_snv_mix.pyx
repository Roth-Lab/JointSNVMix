DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9
DEF NUM_BASES = 2

cdef class JointSnvMixClassifier(Classifier):
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return JointSnvMixClassifierRefIterator(
                                               ref,
                                               self._counter.iter_ref(ref),
                                               **kwargs
                                               )
             
cdef class JointSnvMixClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryBaseCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
    
        self._init_params(**kwargs)
        
    def _init_params(self, **kwargs):
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        
        pi = kwargs.get('pi', (1e6, 1e2, 1e2, 1e2, 1e4, 1e2, 1e1, 1e1, 1e4))
        
        for i in range(NUM_GENOTYPES):
            self._log_mu_N[i][0] = log(< double > mu_N[i])
            self._log_mu_N[i][1] = log(< double > (1 - mu_N[i]))
            
            self._log_mu_T[i][0] = log(< double > mu_T[i])
            self._log_mu_T[i][1] = log(< double > (1 - mu_T[i]))
        
        nc = sum(pi)
        for i in range(NUM_JOINT_GENOTYPES):            
            self._log_pi[i] = log(< double > pi[i] / nc)


    cdef tuple _get_labels(self):
        cdef int  g
        cdef int normal_counts[2], tumour_counts[2]
        cdef double x
        cdef double normal_log_likelihood[NUM_GENOTYPES], tumour_log_likelihood[NUM_GENOTYPES]
        cdef double * normal_probabilities, * tumour_probabilities
        cdef double * joint_probabilities 
        cdef JointBinaryCounterRow row
        cdef tuple labels
        
        row = self._iter._current_row
        
        normal_counts[0] = row._normal_counts.A
        normal_counts[1] = row._normal_counts.B
        
        tumour_counts[0] = row._tumour_counts.A
        tumour_counts[1] = row._tumour_counts.B
        
        for g in range(NUM_GENOTYPES):      
            normal_log_likelihood[g] = multinomial_log_likelihood(normal_counts, self._log_mu_N[g], NUM_BASES)
    
            tumour_log_likelihood[g] = multinomial_log_likelihood(tumour_counts, self._log_mu_T[g], NUM_BASES)

        joint_probabilities = get_joint_posterior(normal_log_likelihood,
                                                  tumour_log_likelihood,
                                                  self._log_pi,
                                                  NUM_GENOTYPES,
                                                  NUM_GENOTYPES)
        
        labels = tuple([x for x in joint_probabilities[:NUM_JOINT_GENOTYPES]])
        
        # Cleanup allocated arrays.
        free(joint_probabilities)
        
        return labels
