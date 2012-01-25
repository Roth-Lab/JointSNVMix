'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

import ConfigParser

cdef class SnvMixTwoModel(MixtureModel):
    def __cinit__(self, BinomialPriors priors, BinomialParameters params):
        self._density = SnvMixTwoDensity(params)
        
        self._ess = SnvMixTwoEss(len(params.mu_N), len(params.mu_T))

    cdef _get_prior_log_likelihood(self):
        '''
        Compute the prior portion of the log likelihood.
        '''        
        ll = 0
        
        for mu_N, mu_N_prior in zip(self.params.mu_N, self.priors.mu_N):
            ll += beta_log_likelihood(mu_N, mu_N_prior['alpha'], mu_N_prior['beta'])

        for mu_T, mu_T_prior in zip(self.params.mu_T, self.priors.mu_T):
            ll += beta_log_likelihood(mu_T, mu_T_prior['alpha'], mu_T_prior['beta'])         
        
        ll += dirichlet_log_likelihood(self.params.pi, self.priors.pi)

        return ll

cdef class SnvMixTwoDensity(Density):
    def __cinit__(self, BinomialParameters params):        
        self._num_normal_genotypes = len(params.mu_N)
        
        self._num_tumour_genotypes = len(params.mu_T)
           
        self._num_joint_genotypes = self._num_normal_genotypes * self._num_tumour_genotypes

        self._init_arrays()

        self.set_params(params)        
        
    def __dealloc__(self):
        free(self._mu_N)
        free(self._mu_T)
        free(self._log_mix_weights)
    
    cdef _init_arrays(self):
        self._mu_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._mu_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        
        self._log_mix_weights = < double *> malloc(sizeof(double) * self._num_joint_genotypes)    

    cdef set_params(self, Parameters params):
        for i, mu_N in enumerate(params.mu_N):
            self._mu_N[i] = mu_N

        for i, mu_T in enumerate(params.mu_T):
            self._mu_T[i] = mu_T

        for i, pi in enumerate(params.pi):
            self._log_mix_weights[i] = log(pi)

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

cdef class SnvMixTwoEss(Ess):    
    def __cinit__(self, int num_normal_genotypes, int num_tumour_genotypes):        
        self._init_arrays()
    
    def __dealloc__(self):
        free(self._a_N)
        free(self._b_N)
        free(self._a_T)
        free(self._b_T)

        free(self._mu_N)
        free(self._mu_T)
    
    cdef _init_arrays(self):
        self._a_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._b_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)

        self._a_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        self._b_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
       
        self._n = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
        
        self._mu_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._mu_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        
    cdef reset(self):
        cdef int i
    
        for i in range(self._num_normal_genotypes):
            self._a_N[i] = 0
            self._b_N[i] = 0
            
        for i in range(self._num_tumour_genotypes):
            self._a_T[i] = 0
            self._b_T[i] = 0            

        for i in range(self._num_joint_genotypes):
            self._n[i] = 0

    cdef set_params(self, Parameters params):
        '''
        Copy Python level parameters into C arrays for fast access.
        '''
        for i, mu_N in enumerate(params.mu_N):
            self._mu_N[i] = mu_N

        for i, mu_T in enumerate(params.mu_T):
            self._mu_T[i] = mu_T

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

    property a_N:
        def __get__(self):
            return [x for x in self._a_N[:self._num_normal_genotypes]]
        
    property a_T:
        def __get__(self):
            return [x for x in self._a_T[:self._num_tumour_genotypes]]
        
    property b_N:
        def __get__(self):
            return [x for x in self._b_N[:self._num_normal_genotypes]]
        
    property b_T:
        def __get__(self):
            return [x for x in self._b_T[:self._num_tumour_genotypes]]  