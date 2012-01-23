'''
Created on 2012-01-22

@author: Andrew Roth
'''
cdef class _JointSnvMixTwoDensity(_JointSnvMixDensity):
    cdef _get_complete_log_likelihood(self, JointBinaryData uncast_data_point, double * ll):
        cdef int g_N, g_T, g_J
        cdef double mu_N, mu_T, log_mix_weight, normal_log_likelihood, tumour_log_likelihood
    
        cdef JointBinaryQualityData data_point = < JointBinaryQualityData > uncast_data_point
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                # Index of joint genotype
                g_J = (self._num_tumour_genotypes * g_N) + g_T
                
                mu_N = self._mu_N[g_N]
                mu_T = self._mu_T[g_T]
                
                log_mix_weight = self._log_mix_weights[g_J]
                                        
                normal_log_likelihood = snv_mix_two_log_likelihood(data_point._q_N,
                                                                   data_point._r_N,
                                                                   data_point._d_N,
                                                                   mu_N)
                
                tumour_log_likelihood = snv_mix_two_log_likelihood(data_point._q_T,
                                                                   data_point._r_T,
                                                                   data_point._d_T,
                                                                   mu_T)
                
                # Combine the mix-weight, normal likelihood and tumour likelihood to obtain class likelihood
                ll[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood


cdef class _JointSnvMixTwoEss(_JointSnvMixEss):
    cdef update(self, JointBinaryData uncast_data_point, double * resp):
        cdef int g_N, g_T, g_J
        cdef double mu_N, mu_T, a_N, a_T, b_N, b_T
        
        cdef JointBinaryQualityData data_point = < JointBinaryQualityData > uncast_data_point
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                # Index of joint genotype
                g_J = (self._num_tumour_genotypes * g_N) + g_T
                
                mu_N = self._mu_N[g_N]
                mu_T = self._mu_T[g_T]
            
                for i in range(data_point._d_N):
                    a_N = snv_mix_two_expected_a(data_point._q_N[i], data_point._r_N[i], mu_N)
                    b_N = snv_mix_two_expected_b(data_point._q_N[i], data_point._r_N[i], mu_N)
                    
                    self._a_N[g_N] += a_N * resp[g_J]
                    self._b_N[g_N] += b_N * resp[g_J]
                    
                for i in range(data_point._d_T):
                    a_T = snv_mix_two_expected_a(data_point._q_T[i], data_point._r_T[i], mu_T)
                    b_T = snv_mix_two_expected_b(data_point._q_T[i], data_point._r_T[i], mu_T)
                    
                    self._a_T[g_T] += a_T * resp[g_J]
                    self._b_T[g_T] += b_T * resp[g_J]                    

                self._n[g_J] += resp[g_J]                