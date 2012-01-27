'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

import ConfigParser

#=======================================================================================================================
# Priors and Parameters
#=======================================================================================================================
cdef class Priors(object):
    def __init__(self, **kwargs):
        default_pi = (2,) * 9
        
        self._pi = tuple(kwargs.get('pi', default_pi))

    def convert_parameter_to_string(self, name, value):
        s = "{0} : ".format(name)        
        s += "\t".join([str(x) for x in value])
        s += "\n"
        
        return s

    def read_from_file(self, file_name):
        '''
        Load priors from a ConfigParser compliant file.
        '''
        pass
    
    property pi:
        def __get__(self):
            return self._pi
    
cdef class Parameters(object):
    def __init__(self, **kwargs):
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        
        self._pi = tuple(kwargs.get('pi', default_pi))
        
        # Normalise pi
        self._pi = tuple([x / sum(self._pi) for x in self._pi])  

    def convert_parameter_to_string(self, name, value):
        s = "{0} : ".format(name)        
        s += "\t".join([str(x) for x in value])
        s += "\n"
        
        return s

    def read_from_file(self, file_name):
        '''
        Read parameters from ConfigParser compliant file.
        '''
        pass
    
    def write_to_file(self, file_name):
        '''
        Write parameters to ConfigParser compliant file.
        '''
        pass
    
    cdef update(self, Ess ess, Priors priors):
        self._update_pi(ess.n, priors.pi)
    
    cdef _update_pi(self, n, prior):
        '''
        Compute the MAP update of the mix-weights in a mixture model with a Dirichlet prior.
        '''        
        pi = []
        
        for n_g, prior_g in zip(n, prior):
            pi.append(n_g + prior_g - 1)
        
        pi = [x / sum(pi) for x in pi]

        self._pi = tuple(pi)      
    
    property pi:
        def __get__(self):
            return self._pi    

cdef class MixtureModel(object):
    #===================================================================================================================
    # Needs to be implemented.
    #===================================================================================================================
    cdef _get_prior_log_likelihood(self):
        '''
        Compute the prior portion of the log likelihood.
        '''        
        pass    

    #===================================================================================================================
    # Implementated
    #===================================================================================================================
    def __cinit__(self, Priors priors, Parameters params):
        self._priors = priors
        self._params = params
        
        self._num_joint_genotypes = len(params.pi)
        
        self._resp = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
    
    def __dealloc__(self):
        free(self._resp)
        
    def predict(self, data_point):
        self._predict(data_point)
        
        return [x for x in self._resp[:self._num_joint_genotypes]]
    
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
        '''
        Update self._ess using data.
        '''
        cdef JointBinaryData data_point
        
        self._ess.reset()
        self._ess.set_params(self._params)
        self._density.set_params(self._params)

        for data_point in data:
            self._density.get_responsibilities(data_point, self._resp)
            self._ess.update(data_point, self._resp)

    cdef _M_step(self):
        '''
        Update self._params.
        '''
        self._params.update(self._ess, self._priors)             
    
    cdef _predict(self, JointBinaryData data_point):
        '''
        C level predict method. After call results are stored in self._resp.
        '''
        self._density.get_responsibilities(data_point, self._resp)
    
    cdef double _get_log_likelihood(self, data):
        '''
        Get log likelihood of data-set.
        '''
        cdef double log_likelihood
        cdef JointBinaryData data_point
        
        log_likelihood = self._get_prior_log_likelihood()
        
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
cdef class Density(object):
    '''
    Base class for density objects. Sub-classing objects need to implement one method, get_responsibilities. This method
    computes the responsibilities for a data-point.
    '''   
    cdef set_params(self, Parameters params):
        '''
        Copy Python level parameters into C arrays for fast access.
        '''
        pass

    cdef _get_complete_log_likelihood(self, JointBinaryData data_point, double * ll):
        '''
        Get the log_likelihood the data point belongs to each class in the model. This will be stored in ll.
        '''
        pass
    
    #===================================================================================================================
    # Implemented
    #===================================================================================================================
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

cdef class Ess(object):
    '''
    Base class for storing and updating expected sufficient statistics (ESS) for JointSnvMix models.
    '''
    def __cinit__(self, int num_normal_genotypes, int num_tumour_genotypes):        
        self._num_normal_genotypes = num_normal_genotypes
        
        self._num_tumour_genotypes = num_tumour_genotypes
        
        self._num_joint_genotypes = num_normal_genotypes * num_tumour_genotypes        

        self._n = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
        
        self.reset()
    
    def __dealloc__(self):
        free(self._n)
    
    #===================================================================================================================
    # Interface
    #===================================================================================================================    
    cdef update(self, JointBinaryData data_point, double * resp):
        '''
        Update the ESS given the data-point and responsibilities.
        '''
        pass
    
        
    cdef reset(self):
        '''
        Reset ESS to 0.
        '''
        pass

    cdef set_params(self, Parameters params):
        '''
        Copy parameters into object.
        '''
        pass
    
    property n:
        def __get__(self):
            return [x for x in self._n[:self._num_joint_genotypes]]
