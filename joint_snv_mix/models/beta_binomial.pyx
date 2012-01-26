'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

import ConfigParser

import joint_snv_mix.constants as constants

cdef class BetaBinomialPriors(Priors):
    def __str__(self):
        return self.convert_parameter_to_string('pi', self.pi)

    def read_from_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        self._pi = tuple([config.getfloat('pi', g) for g in constants.joint_genotypes])

#---------------------------------------------------------------------------------------------------------------------- 
cdef class BetaBinomialParameters(Parameters):
    def __init__(self, **kwargs):
        Parameters.__init__(self, **kwargs)
        
        default_alpha = (1000, 500, 1)        
        self._alpha_N = tuple(kwargs.get('alpha_N', default_alpha))
        self._alpha_T = tuple(kwargs.get('alpha_T', default_alpha))
        
        default_beta = (1, 500, 1000)
        self._beta_N = tuple(kwargs.get('beta_N', default_beta))
        self._beta_T = tuple(kwargs.get('beta_T', default_beta))        
        
    def __str__(self):
        s = ''
        s += self.convert_parameter_to_string('alpha_N', self.alpha_N)
        s += self.convert_parameter_to_string('alpha_T', self.alpha_T)
        s += self.convert_parameter_to_string('beta_N', self.beta_N)
        s += self.convert_parameter_to_string('beta_T', self.beta_T)        
        s += self.convert_parameter_to_string('pi', self.pi)
        
        return s
        
    def read_from_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        self._alpha_N = tuple([config.getfloat('alpha_N', g) for g in constants.genotypes])
        self._alpha_T = tuple([config.getfloat('alpha_T', g) for g in constants.genotypes])
        
        self._beta_N = tuple([config.getfloat('beta_N', g) for g in constants.genotypes])
        self._beta_T = tuple([config.getfloat('beta_T', g) for g in constants.genotypes])

        pi = [config.getfloat('pi', g) for g in constants.joint_genotypes]

        # Normalise pi
        self._pi = tuple([x / sum(pi) for x in pi])

    def write_to_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        
        # Add sections
        config.add_section('pi')
        
        config.add_section('alpha_N')
        config.add_section('alpha_T')
        
        config.add_section('beta_N')
        config.add_section('beta_T')
        
        # Write parameters
        for i, g in enumerate(constants.genotypes):
            config.set('alpha_N', g, str(self.alpha_N[i]))
            config.set('alpha_T', g, str(self.alpha_T[i]))
            
            config.set('beta_N', g, str(self.beta_N[i]))
            config.set('beta_T', g, str(self.beta_T[i]))
                
        for i, g in enumerate(constants.joint_genotypes):
            config.set('pi', g, str(self.pi[i]))
        
        fh = open(file_name, 'w')
        config.write(fh)
        fh.close()

    property alpha_N:
        def __get__(self):
            return self._alpha_N
    
    property alpha_T:
        def __get__(self):
            return self._alpha_T
        
    property beta_N:
        def __get__(self):
            return self._beta_N
    
    property beta_T:
        def __get__(self):
            return self._beta_T        

#=======================================================================================================================
# Model
#=======================================================================================================================
cdef class BetaBinomialModel(MixtureModel):
    def __cinit__(self, BetaBinomialPriors priors, BetaBinomialParameters params):
        self._density = BetaBinomialDensity(params)
        
        self._ess = BetaBinomialEss(len(params.alpha_N), len(params.alpha_T))            
    
    cdef _get_prior_log_likelihood(self):
        ll = 0

        ll += dirichlet_log_likelihood(self.params.pi, self.priors.pi)

        return ll

cdef class BetaBinomialDensity(Density):
    def __cinit__(self, BetaBinomialParameters params):        
        self._num_normal_genotypes = len(params._alpha_N)
        
        self._num_tumour_genotypes = len(params._alpha_T)
           
        self._num_joint_genotypes = self._num_normal_genotypes * self._num_tumour_genotypes

        self._init_arrays()

        self.set_params(params)        
        
    def __dealloc__(self):
        free(self._alpha_N)
        free(self._alpha_T)
        
        free(self._beta_N)
        free(self._beta_T)
        
        free(self._log_mix_weights)

    cdef _init_arrays(self):
        self._alpha_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._alpha_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        
        self._beta_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._beta_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
    
        self._log_mix_weights = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
    
    cdef set_params(self, Parameters params):        
        for i in range(self._num_normal_genotypes):
            self._alpha_N[i] = params.alpha_N[i]
            self._beta_N[i] = params.beta_N[i]

        for i in range(self._num_tumour_genotypes):
            self._alpha_T[i] = params.alpha_T[i]
            self._beta_T[i] = params.beta_T[i]

        for i, pi in enumerate(params.pi):
            self._log_mix_weights[i] = log(pi)
            
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
                
                normal_log_likelihood = beta_binomial_log_likelihood(data_point._a_N, data_point._b_N, alpha_N, beta_N)
                tumour_log_likelihood = beta_binomial_log_likelihood(data_point._a_T, data_point._b_T, alpha_T, beta_T)
                
                # Combine the mix-weight, normal likelihood and tumour likelihood to obtain class likelihood
                ll[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood

cdef class BetaBinomialEss(Ess):
    cdef reset(self):
        cdef int i

        for i in range(self._num_joint_genotypes):
            self._n[i] = 0
    
    cdef set_params(self, Parameters params):
        pass    
    
    cdef update(self, JointBinaryData data_point, double * resp):
        cdef int g_N, g_T, g_J
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                g_J = (self._num_tumour_genotypes * g_N) + g_T

                self._n[g_J] += resp[g_J]

