'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

import ConfigParser

import joint_snv_mix.constants as constants

#=======================================================================================================================
# Priors and Parameters
#=======================================================================================================================
cdef class BinomialPriors(Priors):
    def __init__(self, **kwargs):
        Priors.__init__(self, **kwargs)
        
        default_mu = (
                      {'alpha' : 1000, 'beta' : 2},
                      {'alpha' : 500, 'beta' : 500},
                      {'alpha' : 2, 'beta' : 1000}
                      )
        
        self._mu_N = tuple(kwargs.get('mu_N', default_mu))
        self._mu_T = tuple(kwargs.get('mu_T', default_mu))
        
    def __str__(self):
        s = ''
                
        s += self.convert_parameter_to_string('mu_N_alpha', [x['alpha'] for x in self.mu_N])
        s += self.convert_parameter_to_string('mu_N_beta', [x['beta'] for x in self.mu_N])
        
        s += self.convert_parameter_to_string('mu_T_alpha', [x['alpha'] for x in self.mu_T])
        s += self.convert_parameter_to_string('mu_T_beta', [x['beta'] for x in self.mu_T])
        
        s += self.convert_parameter_to_string('pi', self.pi)
        
        return s

    def read_from_file(self, file_name):       
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        self._mu_N = ({},) * len(constants.genotypes)
        self._mu_T = ({},) * len(constants.genotypes)
        
        for i, g in enumerate(constants.genotypes):
            self._mu_N[i]['alpha'] = config.getfloat('mu_N_alpha', g)
            self._mu_N[i]['beta'] = config.getfloat('mu_N_beta', g)
            
            self._mu_T[i]['alpha'] = config.getfloat('mu_T_alpha', g)
            self._mu_T[i]['beta'] = config.getfloat('mu_T_beta', g)
                
        self._pi = tuple([config.getfloat('pi', g) for g in constants.joint_genotypes])
        
    property mu_N:
        def __get__(self):
            return self._mu_N
    
    property mu_T:
        def __get__(self):
            return self._mu_T

#---------------------------------------------------------------------------------------------------------------------- 
cdef class BinomialParameters(Parameters):
    def __init__(self, **kwargs):
        Parameters.__init__(self, **kwargs)
        
        default_mu = (0.99, 0.5, 0.01)
        
        self._mu_N = tuple(kwargs.get('mu_N', default_mu))
        self._mu_T = tuple(kwargs.get('mu_T', default_mu))
        
    def __str__(self):
        s = ''        
        s += self.convert_parameter_to_string('mu_N', self.mu_N)
        s += self.convert_parameter_to_string('mu_T', self.mu_T)
        s += self.convert_parameter_to_string('pi', self.pi)
        
        return s
    
    def write_to_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        
        config.add_section('pi')
        config.add_section('mu_N')
        config.add_section('mu_T')
        
        for i, g in enumerate(constants.genotypes):
            config.set('mu_N', g, str(self.mu_N[i]))
            config.set('mu_T', g, str(self.mu_T[i]))
            
        for i, g in enumerate(constants.joint_genotypes):
            config.set('pi', g, str(self.pi[i]))
        
        fh = open(file_name, 'w')
        config.write(fh)
        fh.close()
        
    def read_from_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        mu_N = [0, ] * len(constants.genotypes)
        mu_T = [0, ] * len(constants.genotypes)
          
        for i, g in enumerate(constants.genotypes):
            mu_N[i] = config.getfloat('mu_N', g)
            mu_T[i] = config.getfloat('mu_T', g)
        
        pi = [config.getfloat('pi', g) for g in constants.joint_genotypes]       
        
        # Normalise pi
        pi = [x / sum(pi) for x in pi]

        # Save as tuples.
        self._pi = tuple(pi)
        self._mu_N = tuple(mu_N)
        self._mu_T = tuple(mu_T)

    cdef update(self, Ess ess, Priors priors):        
        self._mu_N = self._get_updated_mu(ess.a_N, ess.b_N, priors.mu_N)                
        self._mu_T = self._get_updated_mu(ess.a_T, ess.b_T, priors.mu_T)
        
        self._update_pi(ess.n, priors.pi)

    cdef _get_updated_mu(self, a, b, prior):
        '''
        Compute MAP update to binomial parameter mu with a beta prior.
        '''
        mu = []
        
        for a_g, b_g, prior_g in zip(a, b, prior):
            alpha = a_g + prior_g['alpha'] - 1
            beta = b_g + prior_g['beta'] - 1
            
            denom = alpha + beta

            mu.append(alpha / denom)
        
        return tuple(mu)

    property mu_N:
        def __get__(self):
            return self._mu_N
    
    property mu_T:
        def __get__(self):
            return self._mu_T
        
#=======================================================================================================================
# Model
#=======================================================================================================================
cdef class BinomialModel(MixtureModel):
    def __cinit__(self, BinomialPriors priors, BinomialParameters params):
        self._density = BinomialDensity(params)
        self._ess = BinomialEss(len(params.mu_N), len(params.mu_T))            

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

cdef class BinomialDensity(Density):
    def __cinit__(self, Parameters params):        
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
        '''
        Copy Python level parameters into C arrays for fast access.
        '''
        for i, mu_N in enumerate(params.mu_N):
            self._mu_N[i] = mu_N

        for i, mu_T in enumerate(params.mu_T):
            self._mu_T[i] = mu_T
        
        # Store the log of the mix-weights to speed up computation.
        for i, pi in enumerate(params.pi):
            self._log_mix_weights[i] = log(pi)
            
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
                
                normal_log_likelihood = binomial_log_likelihood(data_point._a_N, data_point._b_N, mu_N)
                tumour_log_likelihood = binomial_log_likelihood(data_point._a_T, data_point._b_T, mu_T)
                
                ll[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood            

cdef class BinomialEss(Ess):
    def __cinit__(self, int num_normal_genotypes, int num_tumour_genotypes):        
        self._init_arrays()
    
    def __dealloc__(self):
        free(self._a_N)
        free(self._b_N)
        free(self._a_T)
        free(self._b_T)
    
    cdef _init_arrays(self):
        self._a_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._b_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)

        self._a_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        self._b_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
       
        self._n = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
        
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
        pass

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
