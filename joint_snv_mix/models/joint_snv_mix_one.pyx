'''
Created on 2012-01-16

@author: Andrew Roth
'''

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixOneEss(JointSnvMixEss):
    '''
    ESS for the JointSNVMix1 model.
    '''
    def __dealloc__(self):
        # Free auxillary data structures.
        free(self._resp)    
    
    cdef _init_auxillary_data_structures(self):
        '''
        Allocate _resp which is used to store the class responsibilities for a data point.
        '''
        # Array for storing class responsibilities for a datapoint.
        self._resp = < double *> malloc(sizeof(double) * self._num_joint_genotypes)    
    
    def update(self, data_point):
        self._compute_responsibilites(data_point)
        self._update_ess(data_point)

    cdef _compute_responsibilities(self, BiAllelicCountDataPoint data_point):
        '''
        Computes the responsibilites and stores them in _resp
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
                        
                normal_log_likelihood = binomial_log_likelihood(data_point._a_N, data_point._b_N, mu_N)
                tumour_log_likelihood = binomial_log_likelihood(data_point._a_T, data_point._b_T, mu_T)
                
                # Combine the mix-weight, normal likelihood and tumour likelihood to obtain class likelihood
                self._resp[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood
        
        # Normalise the class log likelihoods in place to get class posteriors
        normalise_log_space(self._resp, self._num_joint_genotypes)
    
    cdef _update_ess(self, BiAllelicCountDataPoint data_point):
        cdef int g_N, g_T, g_J
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                g_J = (self._tumour_genotypes * g_N) + g_T
            
                self._a_N[g_N] += data_point._a_N * self._resp[g_J]
                self._b_N[g_N] += data_point._b_N * self._resp[g_J]
                
                self._a_T[g_T] += data_point._a_T * self._resp[g_J]
                self._b_T[g_T] += data_point._b_T * self._resp[g_J]
            
                self._n += self._resp[g_J]