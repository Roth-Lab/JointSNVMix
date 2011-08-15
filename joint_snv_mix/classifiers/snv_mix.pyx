#cython: cdivision=True

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9
DEF NUM_BASES = 2

cdef class SnvMixOneClassifier(Classifier):
    def __init__(self, **kwargs):
        self._init_params(**kwargs)
        
    def _init_params(self, **kwargs):
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        for i in range(NUM_GENOTYPES):
            self._log_mu_N[i][0] = log(mu_N[i])
            self._log_mu_N[i][1] = log(1 - mu_N[i])
            
            self._log_mu_T[i][0] = log(mu_T[i])
            self._log_mu_T[i][1] = log(1 - mu_T[i])
            
            self._log_pi_N[i] = log(pi_N[i])            
            self._log_pi_T[i] = log(pi_T[i])

    cdef double * _get_labels(self, PairedSampleBinomialCounterRow row):
        cdef int  g
        cdef int normal_counts[2], tumour_counts[2]
        cdef double x
        cdef double normal_log_likelihood[NUM_GENOTYPES], tumour_log_likelihood[NUM_GENOTYPES]
        cdef double * normal_probabilities, * tumour_probabilities
        cdef double * joint_probabilities
        
        normal_counts[0] = (< JointBinaryCounterRow > row)._normal_counts.A
        normal_counts[1] = (< JointBinaryCounterRow > row)._normal_counts.B
        
        tumour_counts[0] = (< JointBinaryCounterRow > row)._tumour_counts.A
        tumour_counts[1] = (< JointBinaryCounterRow > row)._tumour_counts.B
        
        for g in range(NUM_GENOTYPES):      
            normal_log_likelihood[g] = multinomial_log_likelihood(normal_counts, self._log_mu_N[g], NUM_BASES)
    
            tumour_log_likelihood[g] = multinomial_log_likelihood(tumour_counts, self._log_mu_T[g], NUM_BASES)
        
        normal_probabilities = get_mixture_posterior(normal_log_likelihood, self._log_pi_N, NUM_GENOTYPES)
        tumour_probabilities = get_mixture_posterior(tumour_log_likelihood, self._log_pi_T, NUM_GENOTYPES)
        
        joint_probabilities = combine_independent_probs(normal_probabilities,
                                                        tumour_probabilities,
                                                        NUM_GENOTYPES,
                                                        NUM_GENOTYPES)
        
        # Cleanup allocated arrays.
        free(normal_probabilities)
        free(tumour_probabilities)
        
        return joint_probabilities
