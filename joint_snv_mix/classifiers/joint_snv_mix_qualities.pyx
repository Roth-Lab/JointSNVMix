DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class JointSnvMix2Classifier(Classifier):
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return JointSnvMix2ClassifierRefIterator(
                                                 ref,
                                                 self._counter.iter_ref(ref),
                                                 **kwargs
                                                 )
             
cdef class JointSnvMix2ClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryQualityCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        self._qual_map = get_phred_to_prob_qual_map(256)
        
        self._init_params(**kwargs)

    def __dealloc__(self):
        free(self._qual_map)

    def _init_params(self, **kwargs):
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        
        for i in range(NUM_GENOTYPES):
            self._mu_N[i] = mu_N[i]            
            self._mu_T[i] = mu_T[i]
        
        pi = kwargs.get('pi', (1e6, 1e2, 1e2, 1e2, 1e4, 1e2, 1e1, 1e1, 1e4))        
        nc = sum(pi)        
        for i in range(NUM_JOINT_GENOTYPES):            
            self._log_pi[i] = log(< double > pi[i] / nc)        

    cdef tuple _get_labels(self):
        cdef double x
        cdef double * normal_log_likelihood, * tumour_log_likelihood
        cdef double * joint_probabilities 
        cdef JointBinaryQualityCounterRow row
        cdef tuple labels
        
        row = self._iter._current_row
        
        normal_log_likelihood = snv_mix_2_log_likelihood(row._normal_data, self._mu_N, self._qual_map, NUM_GENOTYPES)
        tumour_log_likelihood = snv_mix_2_log_likelihood(row._tumour_data, self._mu_T, self._qual_map, NUM_GENOTYPES)
        
        joint_probabilities = get_joint_posterior(normal_log_likelihood,
                                                  tumour_log_likelihood,
                                                  self._log_pi,
                                                  NUM_GENOTYPES,
                                                  NUM_GENOTYPES)
        
        labels = tuple([x for x in joint_probabilities[:NUM_JOINT_GENOTYPES]])
        
        # Cleanup allocated arrays.
        free(normal_log_likelihood)
        free(tumour_log_likelihood)
        free(joint_probabilities)
        
        return labels
