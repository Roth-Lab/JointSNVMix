'''
Created on 2012-01-16

@author: Andrew Roth
'''
from __future__ import division

import ConfigParser

#=======================================================================================================================
# Priors and Parameters
#=======================================================================================================================
class JointSnvMixPriors(object):
    def __init__(self, mu_N=None, mu_T=None, pi=None):
        default_mu = (
                      {'alpha' : 100, 'beta' : 2},
                      {'alpha' : 50, 'beta' : 50},
                      {'alpha' : 2, 'beta' : 100}
                      )
        
        default_pi = (2,) * 9
        
        if mu_N is None:
            self.mu_N = default_mu
        
        if mu_T is None:
            self.mu_T = default_mu
        
        if pi is None:
            self.pi = default_pi

    def __str__(self):
        s = "mu_N_alpha : "        
        s += "\t".join([str(x['alpha']) for x in self.mu_N])
        s += "\n"

        s += "mu_N_beta : "        
        s += "\t".join([str(x['beta']) for x in self.mu_N])
        s += "\n"
        
        s += "mu_T_alpha : "        
        s += "\t".join([str(x['alpha']) for x in self.mu_T])
        s += "\n"

        s += "mu_T_beta : "        
        s += "\t".join([str(x['beta']) for x in self.mu_T])
        s += "\n"                   
        
        s += "pi : "
        s += "\t".join([str(x) for x in self.pi])
        s += "\n"
        
        return s

    def read_from_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        joint_genotypes = []
        
        for g_N in genotypes:
            for g_T in genotypes:
                joint_genotypes.append("_".join((g_N, g_T)))
        
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        self.mu_N = []
        self.mu_T = []
        self.pi = []
        
        for g in range(genotypes):
            mu_N = {}
            mu_N['alpha'] = float(config.get('mu_N_alpha', g))
            mu_N['beta'] = float(config.get('mu_N_beta', g)) 
            
            self.mu_N.append(mu_N)
            
            mu_T = {}
            mu_T['alpha'] = float(config.get('mu_T_alpha', g))
            mu_T['beta'] = float(config.get('mu_T_beta', g)) 
            
            self.mu_T.append(mu_T)
        
        for g in range(joint_genotypes):
            pi = float(config.get('pi', g))
            
            self.pi.append(pi)
        
        self.mu_N = tuple(self.mu_N)
        self.mu_T = tuple(self.mu_T)
        
        # Normalise pi
        self.pi = tuple(self.pi)

#---------------------------------------------------------------------------------------------------------------------- 
class JointSnvMixParameters(object):
    def __init__(self, mu_N=None, mu_T=None, pi=None):
        default_mu = (0.99, 0.5, 0.01)
        
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        
        if mu_N is None:
            self.mu_N = default_mu
        else:
            self.mu_N = tuple(mu_N)
            
        if mu_T is None:
            self.mu_T = default_mu
        else:
            self.mu_T = tuple(mu_T)
        
        if pi is None:
            self.pi = default_pi
        else:
            self.pi = pi
        
        # Normalise pi
        self.pi = tuple([x / sum(self.pi) for x in self.pi])         
        
    def __str__(self):
        s = "mu_N : "
        s += "\t".join([str(x) for x in self.mu_N])
        s += "\n"
        
        s += "mu_T : "
        s += "\t".join([str(x) for x in self.mu_T])
        s += "\n"

        s += "pi : "
        s += "\t".join([str(x) for x in self._pi])
        s += "\n"
        
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
        
        self.mu_N = []
        self.mu_T = []
        self.pi = []
        
        for g in range(genotypes):
            mu_N = config.getfloat('mu_N', g)            
            self.mu_N.append(mu_N)
            
            mu_T = config.getfloat('mu_T', g)
            self.mu_T.append(mu_T)
        
        for g in range(joint_genotypes):
            pi = config.getgetfloat('pi', g)
            self.pi.append(pi)
        
        self.mu_N = tuple(self.mu_N)
        self.mu_T = tuple(self.mu_T)
        
        # Normalise pi
        self.pi = tuple([x / sum(self.pi) for x in self.pi])

#=======================================================================================================================
# Model
#=======================================================================================================================
class JointSnvMixModel(object):
    cdef JointSnvMixPriors priors
    cdef JointSnvMixParameters params
    
    cdef _JointSnvMixDensity _density
    cdef _JointSnvMixEss _ess
    
    cdef int _num_genotypes
    cdef double * _resp
    
    def __cinit__(self, priors, params, model="jsm1"):
        self.priors = priors
        self.params = params
        
        if model == "jsm1":
            self._density = _JointSnvMixOneDensity(params)
            self._ess = _JointSnvMixOneEss(len(params.mu_N), len(params.mu_T))
        elif model == "jsm2":
            self._density = _JointSnvMixTwoDensity(params)
            self._ess = _JointSnvMixTwoEss(len(params.mu_N), len(params.mu_T))
        else:
            raise Exception("{0} not a recongnised model. Options are jsm1, jsm2.".format(model))
        
        self._num_joint_genotypes = len(params.mu_N) * len(params.mu_T)
        self._resp = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
    
    def __dealloc__(self):
        free(self._resp)

    def predict(self, data_point):
        self._density.get_responsibilities(data_point, self._resp)
        
        return [x for x in self._resp[:self._num_joint_genotypes]]
    
    def fit(self, data, max_iters=1000, tolerance=1e-6, verbose=False):
        '''
        Fit the model using the EM algorithm.
        '''        
        iters = 0
        ll = [float('-inf')]
        converged = False
        
        while not converged:            
            ess = self._E_step(data)
            self._M_step(data, ess)
            
            ll_iter = self._get_log_likelihood(data)
            
            ll.append(ll_iter)
            
            ll_diff = ll[-1] - ll[-2]
            
            iters += 1
            
            if verbose:
                print "#" * 20
                print iters, ll[-1]
                print self.params
            
            if ll_diff < 0:
                print self.params
                print ll[-1], ll[-2]
                raise Exception('Lower bound decreased.')
            elif ll_diff < tolerance:
                print "Converged"
                converged = True
            elif iters >= max_iters:
                print "Maximum number of iterations exceeded exiting."
                converged = True
            else:
                converged = False
    
    cdef _E_step(self, data):    
        self._ess.reset()
        self._density.set_params(self.params)

        for data_point in data:
            self._density.get_responsibilities(data_point, self._resp)
            self._ess.update(data_point, self._resp)

    cdef _M_step(self):  
        self.params.mu_N = self._get_updated_mu(self._ess.a_N, self._ess.b_N, self.priors.mu_N)
        self.params.mu_T = self._get_updated_mu(self._ess.a_T, self._ess.b_T, self.priors.mu_T)
        
        self.params.pi = self._get_updated_pi(self._ess.n, self.priors.pi)

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
        
        return mu
            
    cdef _get_updated_pi(self, n, prior):
        '''
        Compute the MAP update of the mix-weights in a mixture model with a Dirichlet prior.
        '''        
        pi = []
        
        for n_g, prior_g in zip(n, prior):
            pi.append(n_g + prior_g - 1)
        
        pi = [x / sum(pi) for x in pi]

        return pi
    
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

#=======================================================================================================================
# Density
#=======================================================================================================================
class _JointSnvMixDensity(object):
    '''
    Base class for density objects. Sub-classing objects need to implement one method, get_responsibilities. This method
    computes the responsibilities for a data-point.
    '''
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _mu_N
    cdef double * _mu_T
    cdef double * _log_mix_weights
    #===================================================================================================================
    # Interface
    #===================================================================================================================
    cdef get_responsibilities(self, JointBinaryCountData data_point, double * resp):
        '''
        Computes the responsibilities of the given data-point. Results are stored in resp.
        '''
        pass
    
    #===================================================================================================================
    # Implementation
    #===================================================================================================================
    def __init__(self, params):        
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
    
    cdef set_params(self, JointSnvMixParams params):
        '''
        Copy Python level parameters into C arrays for fast access.
        '''
        for i, mu_N in enumerate(params.mu_N):
            self._mu_N[i] = mu_N

        for i, mu_T in enumerate(params.mu_T):
            self._mu_N[i] = mu_T
        
        # Store the log of the mix-weights to speed up computation.
        for i, pi in enumerate(params.pi):
            self._log_mix_weights[i] = log(pi)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class _JointSnvMixOneDensity(JointSnvMixDensity):
    cdef get_responsibilities(self, JointBinaryCountData data_point, double * resp):
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
                resp[g_J] = log_mix_weight + normal_log_likelihood + tumour_log_likelihood
        
        # Normalise the class log likelihoods in place to get class posteriors
        normalise_log_space(resp, self._num_joint_genotypes)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class _JointSnvMixTwoDensity(JointSnvMixDensity):
    cdef get_responsibilities(self, JointBinaryCountData data_point, double * resp):
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

#=======================================================================================================================
# Ess
#=======================================================================================================================
cdef class _JointSnvMixEss(object):
    '''
    Base class for storing and updating expected sufficient statistics (ESS) for JointSnvMix models using Bernoulli or
    Binomial distributions.
    '''
    cdef int _num_normal_genotypes
    cdef int _num_tumour_genotypes
    cdef int _num_joint_genotypes
    
    cdef double * _a_N
    cdef double * _b_N
    cdef double * _a_T
    cdef double * _b_T
    cdef double * _n
        
    #===================================================================================================================
    # Interface
    #===================================================================================================================
    cdef update(self, JointBinaryCountData data_point, double * resp):
        '''
        Update the ESS given the data-point and responsibilities.
        '''
        pass
    
    #===================================================================================================================
    # Implementation
    #===================================================================================================================
    def __init__(self, int num_normal_genotypes, int num_tumour_genotypes):        
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

    cdef _init_arrays(self):
        '''
        Allocate arrays for sufficient statistics and initialise to 0.        
        '''
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
                    
#---------------------------------------------------------------------------------------------------------------------- 
class JointSnvMixOneEss(JointSnvMixEss):
    cdef update(self, JointBinaryCountData data_point, double * resp):
        cdef int g_N, g_T, g_J
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                g_J = (self._tumour_genotypes * g_N) + g_T
            
                self.a_N[g_N] += data_point._a_N * resp[g_J]
                self.b_N[g_N] += data_point._b_N * resp[g_J]
                
                self.a_T[g_T] += data_point._a_T * resp[g_J]
                self.b_T[g_T] += data_point._b_T * resp[g_J]
            
                self.n += resp[g_J]      

#---------------------------------------------------------------------------------------------------------------------- 
class JointSnvMixTwoEss(JointSnvMixEss):
    cdef update(self, JointBinaryCountData data_point, double * resp):
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
                    
                    self.a_N[g_N] += a_N * self._resp[g_J]
                    self.b_N[g_N] += b_N * self._resp[g_J]
                    
                for i in range(data_point._num_tumour_reads):
                    a_T = snv_mix_two_expected_a(data_point._q_T[i], data_point._r_T[i], mu_T)
                    b_T = snv_mix_two_expected_b(data_point._q_T[i], data_point._r_T[i], mu_T)
                    
                    self.a_T[g_T] += a_T * self._resp[g_J]
                    self.b_T[g_T] += b_T * self._resp[g_J]                    

                self.n += self._resp[g_J] 
