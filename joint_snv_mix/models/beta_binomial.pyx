'''
Created on 2012-01-23

@author: Andrew Roth
'''
#=======================================================================================================================
# Density
#=======================================================================================================================
cdef class _BetaBinomialDensity(Density):    
    cdef _get_complete_log_likelihood(self, JointBinaryData uncast_data_point, double * ll):        
        cdef int g_N, g_T, g_J
        cdef double alpha_N, alpha_T, beta_N, beta_T, log_mix_weight, normal_log_likelihood, tumour_log_likelihood
        
        cdef JointBinaryCountData data_point = < JointBinaryCountData > uncast_data_point
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                # Index of joint genotype
                g_J = (self._num_tumour_genotypes * g_N) + g_T
                
                alpha_N = self._alpha_N[g_N]
                alpha_T = self._alpha_T[g_T]
                
                beta_N = self._beta_N[g_N]
                beta_T = self._beta_T[g_T]
                
                log_mix_weight = self._log_mix_weights[g_J]
                
                normal_log_likelihood = beta_binomial_log_likelihood(data_point._a_N,
                                                                     data_point._b_N,
                                                                     alpha_N,
                                                                     beta_N)
                
                tumour_log_likelihood = beta_binomial_log_likelihood(data_point._a_T,
                                                                     data_point._b_T,
                                                                     alpha_T,
                                                                     beta_T)
                
                # Combine the mix-weight, normal likelihood and tumour likelihood to obtain class likelihood
                ll[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood

#=======================================================================================================================
# ESS
#=======================================================================================================================
cdef class BetaBinomialEss(Ess):
    cdef update(self, JointBinaryData data_point, double * resp):
        cdef int g_N, g_T, g_J
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                g_J = (self._num_tumour_genotypes * g_N) + g_T
                            
                self._n[g_J] += resp[g_J]

#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class BetaBinomialParameters(ModelParameters):
    def __str__(self):
        return self.param_to_string('pi', self.pi)
    
    def write_to_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        joint_genotypes = []
        
        for g_N in genotypes:
            for g_T in genotypes:
                joint_genotypes.append("_".join((g_N, g_T)))
        
        config = ConfigParser.SafeConfigParser()
        
        config.add_section('pi')
        
        for g_J, pi in zip(joint_genotypes, self.pi):
            config.set('pi', g_J, "{0:.10f}".format(pi))
        
        fh = open(file_name, 'w')
        config.write(fh)
        fh.close()
        
    def read_from_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        joint_genotypes = []
        
        for g_N in genotypes:
            for g_T in genotypes:
                joint_genotypes.append("_".join((g_N, g_T)))
        
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)

        pi = []
               
        for g in joint_genotypes:
            pi_g = config.getfloat('pi', g)
            pi.append(pi_g)
       
        self._pi = tuple([x / sum(pi) for x in pi])

#=======================================================================================================================
# Priors
#=======================================================================================================================
cdef class BetaBinomialHyperParameters(ModelHyperParameters):
    def __str__(self):
        return self.param_to_string('pi', self.pi)

    def read_from_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        joint_genotypes = []
        
        for g_N in genotypes:
            for g_T in genotypes:
                joint_genotypes.append("_".join((g_N, g_T)))
        
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        pi = []
        
        for g in joint_genotypes:
            pi_g = float(config.get('pi', g))
            
            pi.append(pi_g)
        
        self._pi = tuple(pi)

    cdef get_log_likelihood(self, params):
        '''
        Compute the prior portion of the log likelihood.
        '''        
        ll = 0  
        
        ll += dirichlet_log_likelihood(params.pi, self.pi)

        return ll
