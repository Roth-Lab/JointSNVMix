'''
Created on 2012-01-16

@author: Andrew Roth
'''
class Parameters(object):
    def __init__(self, pi=None, mu_N=None, mu_T=None, num_genotypes=9):
        if pi is None:
            self.pi = tuple([1 / num_genotypes] * num_genotypes) 
        else:
            self.pi = pi
        
        default_mu = [0.999, 0.5, 0.001]
        
        if mu_N is None:
            self.mu_N = default_mu
        else:
            self.mu_N = mu_N
            
        if mu_T is None:
            self.mu_T = default_mu
        else:
            self.mu_T = mu_T

#=======================================================================================================================
# Model
#=======================================================================================================================
class JointSnvMix(object):
    def __init__(self, priors, params):
        self.priors = priors
        self.params = params
        
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
    
    def _E_step(self, data):
        ess = JointSnvMixEss(self.params)
        
        log_mix_weights = [log(pi) for pi in self.params.pi]
        
        for data_point in data:
            ess.update(data_point)

    def _M_step(self, ess):
        self.params.mu_N = self._get_updated_mu(ess.a_N, ess.b_N, self.priors.mu_N)
        self.params.mu_T = self._get_updated_mu(ess.a_T, ess.b_T, self.priors.mu_T)
        
        self.params.pi = self._get_updated_pi(ess.n, self.priors.pi)

    def _get_updated_mu(self, a, b, prior):
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
            
    def _get_updated_pi(self, n, prior):
        '''
        Compute the MAP update of the mix-weights in a mixture model with a Dirichlet prior.
        '''        
        pi = []
        
        for n_g, prior_g in zip(n, prior):
            pi.append(n_g + prior_g - 1)
        
        pi = [x / sum(pi) for x in pi]

        return pi
    
    def _get_prior_log_likelihood(self):
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
# ESS
#=======================================================================================================================
class JointSnvMixEss(object):
    def __init__(self, params):        
        self._num_normal_genotypes = len(params.mu_N)
        
        self._num_tumour_genotypes = len(params.mu_T)
        
        self._num_joint_genotypes = self._num_normal_genotypes * self._num_tumour_genotypes

        self._copy_params(params)

        self._init_ess()
        
        self._init_auxillary_data_structures()
        
    def __dealloc__(self):
        '''
        Free allocated memory. By Cython convention any subclass will call this in addition to its own __dealloc__
        method.
        '''
        # Free parameter arrays
        free(self._mu_N)
        free(self._mu_T)
        free(self._log_mix_weights)
        
        # Free ess arrays
        free(self._a_N)
        free(self._b_N)
        free(self._a_T)
        free(self._b_T)
    
    cdef _copy_params(self, params):
        '''
        Copy Python level parameters into C arrays for fast access.
        '''
        self._mu_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._mu_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        
        self._log_mix_weights = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
        
        for i, mu_N in enumerate(params.mu_N):
            self._mu_N[i] = mu_N

        for i, mu_T in enumerate(params.mu_T):
            self._mu_N[i] = mu_T
        
        # Store the log of the mix-weights to speed up computation.
        for i, pi in enumerate(params.pi):
            self._log_mix_weights[i] = log(pi)
        
    cdef _init_ess(self):
        '''
        Allocate arrays for sufficient statistics and initialise to 0.        
        '''
        self._a_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        self._b_N = < double *> malloc(sizeof(double) * self._num_normal_genotypes)
        
        for i in range(self._num_normal_genotypes):
            self._a_N[i] = 0
            self._b_N[i] = 0
            
        self._a_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        self._b_T = < double *> malloc(sizeof(double) * self._num_tumour_genotypes)
        
        for i in range(self._num_tumour_genotypes):
            self._a_T[i] = 0
            self._b_T[i] = 0            
       
        self._n = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
        
        for i in range(self._num_joint_genotypes):
            self._n[i] = 0

    cdef _init_auxillary_data_structures(self):
        '''
        Allocate memory for any other data structures which are recurrently used by a model. Any variables allocated
        here should be deallocated in the __dealloc__ method of the subclass.
        '''
        pass
        
    def update(self, data_point):
        '''
        Update the ESS given the data point.
        '''
        pass

cdef class BinomialEss(JointSnvMixEss):
    '''
    ESS for the JointSNVMix1 model.
    '''
    def __dealloc__(self):
        # Free auxillary data structures.
        free(self._resp)    
    
    cdef _init_auxillary_data_structures(self):
        '''
        Allocate _resp which is used to store the class responsibilites for a data point.
        '''
        # Array for storing class responsibilities for a datapoint.
        self._resp = < double *> malloc(sizeof(double) * self._num_joint_genotypes)    
    
    def update(self, data_point):
        self._compute_responsibilites(data_point)
        self._update_ess(data_point)

    cdef _compute_responsibilities(self, DataPoint data_point):
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
    
    cdef _update_ess(self, DataPoint data_point):
        cdef int g_N, g_T, g_J
    
        for g_N in range(self._num_normal_genotypes):            
            for g_T in range(self._num_tumour_genotypes):
                g_J = (self._tumour_genotypes * g_N) + g_T
            
                self._a_N[g_N] += data_point._a_N * self._resp[g_J]
                self._b_N[g_N] += data_point._b_N * self._resp[g_J]
                
                self._a_T[g_T] += data_point._a_T * self._resp[g_J]
                self._b_T[g_T] += data_point._b_T * self._resp[g_J]
            
                self._n += self._resp[g_J]

cdef class JointSnvMixTwoEss(ESS):
    '''
    ESS for the JointSNVMix2 model.
    '''
    cdef _init_auxillary_data_structures(self):
        '''
        Allocate _resp which is used to store the class responsibilites for a data point.
        '''
        # Array for storing class responsibilities for a datapoint.
        self._resp = < double *> malloc(sizeof(double) * self._num_joint_genotypes)
    
    cdef _compute_responsibilities(self, DataPoint data_point):
        '''
        Computes the posteriors and stores them in _resp
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
        
    cdef _update_ess(self, DataPoint data_point):
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
                    
                    self._a_N[g_N] += a_N * self._resp[g_J]
                    self._b_N[g_N] += b_N * self._resp[g_J]
                    
                for i in range(data_point._num_tumour_reads):
                    a_T = snv_mix_two_expected_a(data_point._q_T[i], data_point._r_T[i], mu_T)
                    b_T = snv_mix_two_expected_b(data_point._q_T[i], data_point._r_T[i], mu_T)
                    
                    self._a_T[g_T] += a_T * self._resp[g_J]
                    self._b_T[g_T] += b_T * self._resp[g_J]                    

                self._n += self._resp[g_J]        
