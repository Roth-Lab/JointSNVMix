'''
Created on 2012-01-16

@author: Andrew Roth
'''
cdef class JointSnvMixTwoEss(ESS):
    '''
    ESS for the JointSNVMix2 model.
    '''
    cdef _init_auxillary_data_structures(self):
        '''
        Allocate _resp which is used to store the class responsibilities for a data point.
        '''
        # Array for storing class responsibilities for a datapoint.
        self._resp = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
    
    cdef _compute_responsibilities(self, BiAllelicQualityDataPoint data_point):
        '''
        Computes the posteriors and stores them in _resp
        '''
        cdef int g_N, g_T, g_J
        cdef double mu_N, mu_T, log_mix_weight, normal_log_likelihood, tumour_log_likelihood
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                # Index of joint genotype
                g_J = (self._tumour_genotypes * g_N) + g_T
                
                mu_N = self._mu_N[g_N]
                mu_T = self._mu_T[g_T]
                
                log_mix_weight = self._log_mix_weights[g_J]
                                        
                normal_log_likelihood = snv_mix_two_log_likelihood(data_point._q_N,
                                                                   data_point._r_N,
                                                                   data_point._num_normal_reads,
                                                                   mu_N)
                
                tumour_log_likelihood = snv_mix_two_log_likelihood(data_point._q_T,
                                                                   data_point._r_T,
                                                                   data_point._num_tumour_reads,
                                                                   mu_T)
                
                # Combine the mix-weight, normal likelihood and tumour likelihood to obtain class likelihood
                self._resp[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood
        
        # Normalise the class log likelihoods in place to get class posteriors
        normalise_log_space(self._resp, self._num_joint_genotypes)
        
    cdef _update_ess(self, BiAllelicQualityDataPoint data_point):
        cdef int g_N, g_T, g_J
        cdef double mu_N, mu_T, a_N, a_T, b_N, b_T
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                # Index of joint genotype
                g_J = (self._tumour_genotypes * g_N) + g_T
                
                mu_N = self._mu_N[g_N]
                mu_T = self._mu_T[g_T]
            
                for i in range(data_point._num_normal_reads):
                    a_N = snv_mix_two_expected_a(data_point._q_N[i], data_point._r_N[i], mu_N)
                    b_N = snv_mix_two_expected_b(data_point._q_N[i], data_point._r_N[i], mu_N)
                    
                    self._a_N[g_N] += a_N * self._resp[g_J]
                    self._b_N[g_N] += b_N * self._resp[g_J]
                    
                for i in range(data_point._num_tumour_reads):
                    a_T = snv_mix_two_expected_a(data_point._q_T[i], data_point._r_T[i], mu_T)
                    b_T = snv_mix_two_expected_b(data_point._q_T[i], data_point._r_T[i], mu_T)
                    
                    self._a_T[g_T] += a_T * self._resp[g_J]
                    self._b_T[g_T] += b_T * self._resp[g_J]                    

                self._n += self._resp[g_J]  