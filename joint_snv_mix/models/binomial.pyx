'''
Created on 2012-01-23

@author: Andrew Roth
'''
import ConfigParser

#=======================================================================================================================
# Density
#=======================================================================================================================
cdef class BinomialDensity(Density):    
    cdef _get_complete_log_likelihood(self, JointBinaryData uncast_data_point, double * ll):        
        cdef int g_N, g_T, g_J
        cdef double mu_N, mu_T, log_mix_weight, normal_log_likelihood, tumour_log_likelihood
        
        cdef JointBinaryCountData data_point = < JointBinaryCountData > uncast_data_point
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                # Index of joint genotype
                g_J = (self._num_tumour_genotypes * g_N) + g_T
                
                mu_N = self._mu_N[g_N]
                mu_T = self._mu_T[g_T]
                
                log_mix_weight = self._log_mix_weights[g_J]
                
                normal_log_likelihood = binomial_log_likelihood(data_point._a_N, data_point._b_N, mu_N)
                tumour_log_likelihood = binomial_log_likelihood(data_point._a_T, data_point._b_T, mu_T)
                
                # Combine the mix-weight, normal likelihood and tumour likelihood to obtain class likelihood
                ll[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood

#=======================================================================================================================
# ESS
#=======================================================================================================================
cdef class BinomialEss(Ess):
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

#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class BinomialParameters(ModelParameters):
    def __init__(self, mu_N=None, mu_T=None, pi=None):
        ModelParameters.__init__(self, pi=pi)
        
        default_mu = (0.99, 0.5, 0.01)
        
        if mu_N is None:
            self._mu_N = default_mu
        else:
            self._mu_N = tuple(mu_N)
            
        if mu_T is None:
            self._mu_T = default_mu
        else:
            self._mu_T = tuple(mu_T)
        
    def __str__(self):
        s = ""
        
        s += self.param_to_string("mu_N", self.mu_N)
        s += self.param_to_string("mu_T", self.mu_T)
        s += self.param_to_string("pi", self.pi)
                
        return s
    
    def write_to_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        joint_genotypes = []
        
        for g_N in genotypes:
            for g_T in genotypes:
                joint_genotypes.append("_".join((g_N, g_T)))
        
        config = ConfigParser.SafeConfigParser()
        
        config.add_section('pi')
        config.add_section('mu_N')
        config.add_section('mu_T')
        
        for g_N, mu_N in zip(genotypes, self.mu_N):
            config.set('mu_N', g_N, "{0:.10f}".format(mu_N))
        
        for g_T, mu_T in zip(genotypes, self.mu_T):
            config.set('mu_T', g_T, "{0:.10f}".format(mu_T))
            
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
        
        mu_N = []
        mu_T = []
        pi = []
        
        for g in genotypes:
            mu_N_g = config.getfloat('mu_N', g)            
            mu_N.append(mu_N_g)
            
            mu_T_g = config.getfloat('mu_T', g)
            mu_T.append(mu_T_g)
        
        for g in joint_genotypes:
            pi_g = config.getfloat('pi', g)
            pi.append(pi_g)
        
        self._mu_N = tuple(mu_N)
        self._mu_T = tuple(mu_T)
        
        # Normalise pi
        self._pi = tuple([x / sum(pi) for x in pi])

    property mu_N:
        def __get__(self):
            return self._mu_N
    
    property mu_T:
        def __get__(self):
            return self._mu_T

#=======================================================================================================================
# Priors
#=======================================================================================================================
cdef class BinomialHyperParameters(ModelHyperParameters):
    def __init__(self, mu_N=None, mu_T=None, pi=None):
        ModelHyperParameters.__init__(self, pi=pi)
        
        default_mu = (
                      {'alpha' : 100, 'beta' : 2},
                      {'alpha' : 50, 'beta' : 50},
                      {'alpha' : 2, 'beta' : 100}
                      )
                
        if mu_N is None:
            self._mu_N = default_mu
        else:
            self._mu_N = tuple(mu_N)
        
        if mu_T is None:
            self._mu_T = default_mu
        else:
            self._mu_T = tuple(mu_T)

    def __str__(self):
        s = ""
        
        s += self.param_to_string("mu_N_alpha", [x['alpha'] for x in self.mu_N])
        s += self.param_to_string("mu_N_beta", [x['beta'] for x in self.mu_N])
        
        s += self.param_to_string("mu_T_alpha", [x['alpha'] for x in self.mu_T])
        s += self.param_to_string("mu_T_beta", [x['beta'] for x in self.mu_T])
        
        s += self.param_to_string("pi", self.pi)
                
        return s        

    def read_from_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        joint_genotypes = []
        
        for g_N in genotypes:
            for g_T in genotypes:
                joint_genotypes.append("_".join((g_N, g_T)))
        
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        mu_N = []
        mu_T = []
        pi = []
        
        for g in genotypes:
            mu_N_g = {}
            mu_N_g['alpha'] = float(config.get('mu_N_alpha', g))
            mu_N_g['beta'] = float(config.get('mu_N_beta', g)) 
            
            mu_N.append(mu_N_g)
            
            mu_T_g = {}
            mu_T_g['alpha'] = float(config.get('mu_T_alpha', g))
            mu_T_g['beta'] = float(config.get('mu_T_beta', g)) 
            
            mu_T.append(mu_T_g)
        
        for g in joint_genotypes:
            pi_g = float(config.get('pi', g))
            
            pi.append(pi_g)
        
        self._mu_N = tuple(mu_N)
        self._mu_T = tuple(mu_T)
        
        self._pi = tuple(pi)

    cdef get_log_likelihood(self, params):
        '''
        Compute the prior portion of the log likelihood.
        '''        
        ll = 0
        
        for mu_N, mu_N_prior in zip(params.mu_N, self.mu_N):
            ll += beta_log_likelihood(mu_N, mu_N_prior['alpha'], mu_N_prior['beta'])

        for mu_T, mu_T_prior in zip(params.mu_T, self.mu_T):
            ll += beta_log_likelihood(mu_T, mu_T_prior['alpha'], mu_T_prior['beta'])         
        
        ll += dirichlet_log_likelihood(params.pi, self.pi)

        return ll
        
    property mu_N:
        def __get__(self):
            return self._mu_N
    
    property mu_T:
        def __get__(self):
            return self._mu_T     
