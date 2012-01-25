'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

import ConfigParser

import joint_snv_mix.constants as constants

cdef class BetaBinomialPriors(object):
    def __init__(self, pi=None):
        default_pi = (2,) * 9

        if pi is None:
            self._pi = default_pi
        else:
            self._pi = tuple(pi)

    def __str__(self):
        s = "pi : "
        s += "\t".join([str(x) for x in self.pi])
        s += "\n"
        
        return s

    def read_from_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        pi = []
        
        for g in constants.joint_genotypes:
            pi_g = float(config.get('pi', g))
            
            pi.append(pi_g)
        
        self._pi = tuple(pi)
        
    property pi:
        def __get__(self):
            return self._pi

#---------------------------------------------------------------------------------------------------------------------- 
cdef class BetaBinomialParameters(object):
    def __init__(self, **kwargs):
        default_alpha = (100, 50, 1)        
        self._alpha_N = tuple(kwargs.get('alpha_N', default_alpha))
        self._alpha_T = tuple(kwargs.get('alpha_T', default_alpha))
        
        default_beta = (1, 50, 100)
        self._beta_N = tuple(kwargs.get('beta_N', default_beta))
        self._beta_T = tuple(kwargs.get('beta_T', default_beta))
        
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        self._pi = tuple(kwargs.get('pi', default_pi))
        
        # Normalise pi
        self._pi = tuple([x / sum(self._pi) for x in self._pi])         
        
    def __str__(self):
        s = "pi : "
        s += "\t".join([str(x) for x in self.pi])
        s += "\n"
        
        return s
    
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
            config.set('alpha_N', g, self._alpha_N[i])
            config.set('alpha_T', g, self._alpha_T[i])
            
            config.set('beta_N', g, self._beta_N[i])
            config.set('beta_T', g, self._beta_T[i])
                
        for i, g in enumerate(constants.joint_genotypes):
            config.set('pi', g, self._pi[i])
        
        fh = open(file_name, 'w')
        config.write(fh)
        fh.close()
        
    def read_from_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        self._alpha_N = (0,) * len(constants.genotypes)
        self._alpha_T = (0,) * len(constants.genotypes)
        
        self._beta_N = (0,) * len(constants.genotypes)
        self._beta_T = (0,) * len(constants.genotypes)

        for i, g in enumerate(constants.genotypes):
            self._alpha_N[i] = config.getfloat('alpha_N', g)
            self._alpha_T[i] = config.getfloat('alpha_T', g)
            
            self._beta_N[i] = config.getfloat('beta_N', g)
            self._beta_T[i] = config.getfloat('beta_T', g)
        
        pi = (0,) * len(constants.joint_genotypes)
        
        for i, g in enumerate(constants.joint_genotypes):
            pi[i] = config.getfloat('pi', g)

        # Normalise pi
        self._pi = tuple([x / sum(pi) for x in pi])

    property pi:
        def __get__(self):
            return self._pi
        
#=======================================================================================================================
# Model
#=======================================================================================================================
cdef class BetaBinomialModel(object):
    def __cinit__(self, BetaBinomialPriors priors, BetaBinomialParameters params):
        self._priors = priors
        self._params = params
        
        self._num_normal_genotypes = len(params._alpha_N)
        self._num_tumour_genotypes = len(params._alpha_T)
        self._num_joint_genotypes = self._num_normal_genotypes * self._num_tumour_genotypes
                
        self._density = _Density(params)
        self._ess = _Ess(self._num_normal_genotypes, self._num_tumour_genotypes)            

        self._resp = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
    
    def __dealloc__(self):
        free(self._resp)

    def predict(self, data_point):
        self._predict(data_point)
        
        return [x for x in self._resp[:self._num_joint_genotypes]]
    
    cdef _predict(self, JointBinaryData data_point):
        self._density.get_responsibilities(data_point, self._resp)
    
    def fit(self, data, max_iters=1000, tolerance=1e-6, verbose=False):
        '''
        Fit the model using the EM algorithm.
        '''
        if verbose:        
            print "Fitting model to {0} data-points".format(len(data))
            print
            print "Priors : "
            print self.priors
            print "Initial parameters : "
            print self.params
        
        iters = 0
        ll = [float('-inf')]
        converged = False
        
        while not converged:            
            self._E_step(data)
            
            ll_iter = self._get_log_likelihood(data)
            
            self._M_step()
            
            ll.append(ll_iter)
            
            ll_diff = ll[-1] - ll[-2]
            relative_ll_diff = (ll_diff) / abs(ll[-2]) 
            
            iters += 1
            
            if verbose:
                print "#" * 80
                print iters, ll[-1], ll_diff, relative_ll_diff
                print self.params
            
            if relative_ll_diff < 0:
                print self.params
                print ll[-1], ll[-2], ll_diff, relative_ll_diff
                raise Exception('Lower bound decreased.')
            elif relative_ll_diff < tolerance:
                print "Converged"
                converged = True
            elif iters >= max_iters:
                print "Maximum number of iterations exceeded exiting."
                converged = True
            else:
                converged = False
    
    cdef _E_step(self, data):
        cdef JointBinaryData data_point
        
        self._ess.reset()
        self._density.set_params(self._params)

        for data_point in data:
            self._density.get_responsibilities(data_point, self._resp)
            self._ess.update(data_point, self._resp)

    cdef _M_step(self):
        self._params._pi = self._get_updated_pi(self._ess.n, self._priors._pi)
            
    cdef _get_updated_pi(self, n, prior):
        '''
        Compute the MAP update of the mix-weights in a mixture model with a Dirichlet prior.
        '''        
        pi = []
        
        for n_g, prior_g in zip(n, prior):
            pi.append(n_g + prior_g - 1)
        
        pi = [x / sum(pi) for x in pi]

        return tuple(pi)
    
    cdef double _get_log_likelihood(self, data):
        cdef double log_liklihood
        cdef JointBinaryData data_point
        
        log_likelihood = self._get_prior_log_likelihood()
        
        for data_point in data:
            log_likelihood += self._density.get_log_likelihood(data_point, self._resp)
        
        return log_likelihood
    
    cdef _get_prior_log_likelihood(self):
        ll = 0

        ll += dirichlet_log_likelihood(self.params.pi, self.priors.pi)

        return ll
        
    property params:
        def __get__(self):
            return self._params
        
    property priors:
        def __get__(self):
            return self._priors

#=======================================================================================================================
# Private helper classes
#=======================================================================================================================
cdef class _Density(object):
    def __cinit__(self, BetaBinomialParameters params):        
        self._num_normal_genotypes = len(params._alpha_N)
        
        self._num_tumour_genotypes = len(params._alpha_T)
           
        self._num_joint_genotypes = self._num_normal_genotypes * self._num_tumour_genotypes

        self._init_arrays()

        self.set_params(params)        
        
    def __dealloc__(self):
        free(self._log_mix_weights)
    
    cdef _init_arrays(self):
        self._alpha_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._alpha_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        
        self._beta_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._beta_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
    
        self._log_mix_weights = < double *> malloc(sizeof(double) * self._num_joint_genotypes)

    cdef get_responsibilities(self, JointBinaryData data_point, double * resp):
        cdef int i
        
        self._get_complete_log_likelihood(data_point, resp)
       
        # Normalise the class log likelihoods in place to get class posteriors
        log_space_normalise(resp, self._num_joint_genotypes)
        
        for i in range(self._num_joint_genotypes):
            resp[i] = exp(resp[i])
        
    cdef double get_log_likelihood(self, JointBinaryData data_point, double * ll):
        self._get_complete_log_likelihood(data_point, ll)
        
        return log_sum_exp(ll, self._num_joint_genotypes)
    
    cdef set_params(self, BetaBinomialParameters params):        
        for i in range(self._num_normal_genotypes):
            self._alpha_N[i] = params._alpha_N[i]
            self._beta_N[i] = params._beta_N[i]
        
        for i in range(self._num_tumour_genotypes):
            self._alpha_T[i] = params._alpha_T[i]
            self._beta_T[i] = params._beta_T[i]

        for i, pi in enumerate(params._pi):
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

cdef class _Ess(object):
    def __cinit__(self, int num_normal_genotypes, int num_tumour_genotypes):        
        self._num_normal_genotypes = num_normal_genotypes
        
        self._num_tumour_genotypes = num_tumour_genotypes
        
        self._num_joint_genotypes = num_normal_genotypes * num_tumour_genotypes        

        self._n = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
        
        self.reset()    
    
    def __dealloc__(self):
        free(self._n)
                
    cdef reset(self):
        cdef int i

        for i in range(self._num_joint_genotypes):
            self._n[i] = 0
    
    cdef update(self, JointBinaryData data_point, double * resp):
        cdef int g_N, g_T, g_J
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                g_J = (self._num_tumour_genotypes * g_N) + g_T

                self._n[g_J] += resp[g_J]

