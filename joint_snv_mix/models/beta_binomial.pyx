'''
Created on 2012-01-22

@author: Andrew Roth
'''
cdef class _BetaBinomialDensity(_JointSnvMixDensity):    
    cdef _get_complete_log_likelihood(self, JointBinaryData uncast_data_point, double * ll):        
        cdef int g_N, g_T, g_J, a, b
        cdef double mu_N, mu_T, log_mix_weight, normal_log_likelihood, tumour_log_likelihood
        
        cdef JointBinaryCountData data_point = < JointBinaryCountData > uncast_data_point
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                # Index of joint genotype
                g_J = (self._num_tumour_genotypes * g_N) + g_T
                
                mu_N = self._mu_N[g_N]
                mu_T = self._mu_T[g_T]
                
                log_mix_weight = self._log_mix_weights[g_J]
                
                normal_log_likelihood = beta_binomial_log_likelihood(data_point._a_N, data_point._b_N, mu_N)
                tumour_log_likelihood = beta_binomial_log_likelihood(data_point._a_T, data_point._b_T, mu_T)
                
                # Combine the mix-weight, normal likelihood and tumour likelihood to obtain class likelihood
                ll[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood
                
cdef class _BetaBinomialEss(_JointSnvMixEss):
    cdef update(self, JointBinaryData data_point, double * resp):
        cdef int g_N, g_T, g_J
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                g_J = (self._num_tumour_genotypes * g_N) + g_T
            
                self._a_N[g_N] += data_point._a_N * resp[g_J]
                self._b_N[g_N] += data_point._b_N * resp[g_J]
                
                self._a_T[g_T] += data_point._a_T * resp[g_J]
                self._b_T[g_T] += data_point._b_T * resp[g_J]
            
                self._n[g_J] += resp[g_J]                   