'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division
        
cdef class JointSnvMixModel(object):
    def __cinit__(self, JointSnvMixPriors priors, JointSnvMixParameters params, model="jsm1"):
        self._priors = priors
        self._params = params
        
        if model == "jsm1":
            self._density = _JointSnvMixOneDensity(params)
            self._ess = _JointSnvMixOneEss(len(params._mu_N), len(params._mu_T))            
        elif model == "jsm2":
            self._density = _JointSnvMixTwoDensity(params)
            self._ess = _JointSnvMixTwoEss(len(params._mu_N), len(params._mu_T))
        else:
            raise Exception("{0} not a recongnised model. Options are jsm1, jsm2.".format(model))
        
        self._num_joint_genotypes = len(params._mu_N) * len(params._mu_T)
        
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
        self._ess.set_params(self._params)
        self._density.set_params(self._params)

        for data_point in data:
            self._density.get_responsibilities(data_point, self._resp)
            self._ess.update(data_point, self._resp)

    cdef _M_step(self):
        self._params._mu_N = self._get_updated_mu(self._ess.a_N, self._ess.b_N, self._priors._mu_N)
        self._params._mu_T = self._get_updated_mu(self._ess.a_T, self._ess.b_T, self._priors._mu_T)
        
        self._params._pi = self._get_updated_pi(self._ess.n, self._priors._pi)

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
        
        log_likelihood = self._prior.get_log_likelihood(self._params)
        
        for data_point in data:
            log_likelihood += self._density.get_log_likelihood(data_point, self._resp)
        
        return log_likelihood
    
    property params:
        def __get__(self):
            return self._params
        
    property priors:
        def __get__(self):
            return self._priors

#=======================================================================================================================
# Density
#=======================================================================================================================
cdef class _JointSnvMixDensity(object):
    '''
    Base class for density objects. Sub-classing objects need to implement one method, get_responsibilities. This method
    computes the responsibilities for a data-point.
    '''
    #===================================================================================================================
    # Interface
    #===================================================================================================================
    cdef _get_complete_log_likelihood(self, JointBinaryData data_point, double * ll):
        '''
        Get the log_likelihood the data point belongs to each class in the model. This will be stored in ll.
        '''
        pass
    
    #===================================================================================================================
    # Implementation
    #===================================================================================================================
    def __cinit__(self, JointSnvMixParameters params):        
        self._num_normal_genotypes = len(params._mu_N)
        
        self._num_tumour_genotypes = len(params._mu_T)
           
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

    cdef get_responsibilities(self, JointBinaryData data_point, double * resp):
        '''
        Computes the responsibilities of the given data-point. Results are stored in resp.
        '''
        cdef int i
        
        self._get_complete_log_likelihood(data_point, resp)
       
        # Normalise the class log likelihoods in place to get class posteriors
        log_space_normalise(resp, self._num_joint_genotypes)
        
        for i in range(self._num_joint_genotypes):
            resp[i] = exp(resp[i])
        
    cdef double get_log_likelihood(self, JointBinaryData data_point, double * ll):
        '''
        Computes the log_likelihood for a single point.
        '''
        self._get_complete_log_likelihood(data_point, ll)
        
        return log_sum_exp(ll, self._num_joint_genotypes)
    
    cdef set_params(self, JointSnvMixParameters params):
        '''
        Copy Python level parameters into C arrays for fast access.
        '''
        for i, mu_N in enumerate(params._mu_N):
            self._mu_N[i] = mu_N

        for i, mu_T in enumerate(params._mu_T):
            self._mu_T[i] = mu_T
        
        # Store the log of the mix-weights to speed up computation.
        for i, pi in enumerate(params._pi):
            self._log_mix_weights[i] = log(pi)

#=======================================================================================================================
# Ess
#=======================================================================================================================
cdef class _JointSnvMixEss(object):
    '''
    Base class for storing and updating expected sufficient statistics (ESS) for JointSnvMix models using Bernoulli or
    Binomial distributions.
    '''
    #===================================================================================================================
    # Interface
    #===================================================================================================================
    cdef update(self, JointBinaryData data_point, double * resp):
        '''
        Update the ESS given the data-point and responsibilities.
        '''
        pass
    
    #===================================================================================================================
    # Implementation
    #===================================================================================================================
    def __cinit__(self, int num_normal_genotypes, int num_tumour_genotypes):        
        self._num_normal_genotypes = num_normal_genotypes
        
        self._num_tumour_genotypes = num_tumour_genotypes
        
        self._num_joint_genotypes = num_normal_genotypes * num_tumour_genotypes        

        self._init_arrays()
        
        self.reset()    
    
    def __dealloc__(self):
        free(self._a_N)
        free(self._b_N)
        free(self._a_T)
        free(self._b_T)
        free(self._n)

        free(self._mu_N)
        free(self._mu_T)
    
    cdef _init_arrays(self):
        '''
        Allocate arrays for sufficient statistics and initialise to 0.        
        '''
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

    cdef set_params(self, JointSnvMixParameters params):
        '''
        Copy Python level parameters into C arrays for fast access.
        '''
        for i, mu_N in enumerate(params._mu_N):
            self._mu_N[i] = mu_N

        for i, mu_T in enumerate(params._mu_T):
            self._mu_T[i] = mu_T         
    
    property a_N:
        def __get__(self):
            return [x for x in self._a_N[:self._num_normal_genotypes]]

    property b_N:
        def __get__(self):
            return [x for x in self._b_N[:self._num_normal_genotypes]]
        
    property a_T:
        def __get__(self):
            return [x for x in self._a_T[:self._num_tumour_genotypes]]
        
    property b_T:
        def __get__(self):
            return [x for x in self._b_T[:self._num_tumour_genotypes]]
    
    property n:
        def __get__(self):
            return [x for x in self._n[:self._num_joint_genotypes]]

#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class Parameters(object):
    def param_to_string(self, param_name, param_value):
        s = "{0} : ".format(param_name)
        s += "\t".join([str(x) for x in param_value])
        s += "\n" 
        
    def read_from_file(self, file_name):
        '''
        Read parameters from ConfigParser compliant file.
        '''
        pass        

    property pi:
        def __get__(self):
            return self._pi
                
cdef class ModelParameters(Parameters):
    def __init__(self, pi=None):
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        
        if pi is None:
            self._pi = default_pi
        else:
            self._pi = tuple(pi)
        
        # Normalise pi
        self._pi = tuple([x / sum(self._pi) for x in self._pi])         

    def write_to_file(self, file_name):
        '''
        Write parameters to ConfigParser compliant file.
        '''
        pass
  
cdef class ModelHyperParameters(Parameters):
    def __init__(self, pi=None):
        default_pi = (2,) * 9
        
        if pi is None:
            self._pi = default_pi
        else:
            self._pi = tuple(pi)
            
    def read_from_file(self, file_name):
        '''
        Load parameters from a file compliant with the ConfigParser module.
        '''        
        pass

    cdef get_log_likelihood(self, params):
        '''
        Compute the prior portion of the log likelihood.
        '''        
        pass           
