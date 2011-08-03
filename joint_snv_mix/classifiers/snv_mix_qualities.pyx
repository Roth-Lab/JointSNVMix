#cython: cdivision=True

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMix2Classifier(Classifier):
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return SnvMix2ClassifierRefIterator(
                                           ref,
                                           self._counter.iter_ref(ref),
                                           **kwargs
                                           )
             
cdef class SnvMix2ClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryQualityCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        self._qual_map = get_phred_to_prob_qual_map(256)
        
        self._init_params(**kwargs)
            
    def __dealloc__(self):
        free(self._qual_map)
    
    def _init_params(self, **kwargs):
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        for i in range(NUM_GENOTYPES):
            self._mu_N[i] = mu_N[i]            
            self._mu_T[i] = mu_T[i]
            
            self._log_pi_N[i] = log(pi_N[i])            
            self._log_pi_T[i] = log(pi_T[i])

    cdef tuple _get_labels(self):
        cdef double x
        cdef double * normal_log_likelihood, * tumour_log_likelihood
        cdef double * normal_probabilities, * tumour_probabilities
        cdef double * joint_probabilities 
        cdef JointBinaryQualityCounterRow row
        cdef tuple labels
        
        row = self._iter._current_row
        
        normal_log_likelihood = snv_mix_2_log_likelihood(row._normal_data, self._mu_N, self._qual_map, NUM_GENOTYPES)
        tumour_log_likelihood = snv_mix_2_log_likelihood(row._tumour_data, self._mu_T, self._qual_map, NUM_GENOTYPES)
        
        normal_probabilities = get_mixture_posterior(normal_log_likelihood, self._log_pi_N, NUM_GENOTYPES)
        tumour_probabilities = get_mixture_posterior(tumour_log_likelihood, self._log_pi_T, NUM_GENOTYPES)
        
        joint_probabilities = combine_independent_probs(normal_probabilities,
                                                        tumour_probabilities,
                                                        NUM_GENOTYPES,
                                                        NUM_GENOTYPES)
        
        labels = tuple([x for x in joint_probabilities[:NUM_JOINT_GENOTYPES]])
        
        # Cleanup allocated arrays.
        free(normal_log_likelihood)
        free(tumour_log_likelihood)
        free(normal_probabilities)
        free(tumour_probabilities)
        free(joint_probabilities)
        
        return labels
        
