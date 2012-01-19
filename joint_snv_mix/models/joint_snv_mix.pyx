'''
Created on 2012-01-16

@author: Andrew Roth
'''
#=======================================================================================================================
# Priors
#=======================================================================================================================
cdef class JointSnvMixPriors(object):
    def __init__(self, **kwargs):
        default_mu = (
                      (100, 2),
                      (50, 50),
                      (2, 100)
                      )
        
        default_pi = (2,) * 9
        
        mu_N = kwargs.get('mu_N', default_mu)
        mu_T = kwargs.get('mu_T', default_mu)
        
        pi = kwargs.get('pi', default_pi)
        
        for g in range(NUM_GENOTYPES):
            self._mu_N[g][0] = mu_N[g][0]
            self._mu_N[g][1] = mu_N[g][1]
            
            self._mu_T[g][0] = mu_T[g][0]
            self._mu_T[g][1] = mu_T[g][1]
        
        for g in range(NUM_JOINT_GENOTYPES):    
            self._pi[g] = pi[g]

    def __str__(self):
        s = "mu_N_alpha : "        
        s += "\t".join([str(x[0]) for x in self.mu_N])
        s += "\n"

        s += "mu_N_beta : "        
        s += "\t".join([str(x[1]) for x in self.mu_N])
        s += "\n"
        
        s += "mu_T_alpha : "        
        s += "\t".join([str(x[0]) for x in self.mu_T])
        s += "\n"

        s += "mu_T_beta : "        
        s += "\t".join([str(x[1]) for x in self.mu_T])
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
        
        for g in range(NUM_GENOTYPES):
            self._mu_N[g][0] = float(config.get('mu_N_alpha', genotypes[g]))
            self._mu_T[g][0] = float(config.get('mu_T_alpha', genotypes[g]))
            
            self._mu_N[g][1] = float(config.get('mu_N_beta', genotypes[g]))
            self._mu_T[g][1] = float(config.get('mu_T_beta', genotypes[g]))
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = float(config.get('pi', joint_genotypes[g]))
    
    property mu_N:
        def __get__(self):
            return self._get_mu_as_tuple(self._mu_N)

    property mu_T:
        def __get__(self):
            return self._get_mu_as_tuple(self._mu_T)
    
    property pi:
        def __get__(self):
            return tuple([x for x in self._pi[:NUM_JOINT_GENOTYPES]])
    
    cdef _get_mu_as_tuple(self, double c_mu[3][2]):
            mu = []
            
            for g in range(NUM_GENOTYPES):
                mu.append((c_mu[g][0], c_mu[g][1]))
            
            return tuple(mu)

#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class JointSnvMixParameters(object):
    def __init__(self, **kwargs):        
        self._priors = kwargs.get('priors', JointSnvMixPriors())
        
        default_mu = (0.99, 0.5, 0.01)
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        
        mu_N = kwargs.get('mu_N', default_mu)
        mu_T = kwargs.get('mu_T', default_mu)
        
        pi = kwargs.get('pi', default_pi)
        
        for g in range(NUM_GENOTYPES):
            self._mu_N[g] = mu_N[g]            
            self._mu_T[g] = mu_T[g]
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = pi[g]
        
        self._normalise_pi()            
        
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
        
        for g in range(NUM_GENOTYPES):
            self._mu_N[g] = float(config.get('mu_N', genotypes[g]))
            self._mu_T[g] = float(config.get('mu_T', genotypes[g]))
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = float(config.get('pi', joint_genotypes[g]))
        
        self._normalise_pi()
        
    property mu_N:
        def __get__(self):
            return tuple([x for x in self._mu_N[:NUM_GENOTYPES]])
    
    property mu_T:
        def __get__(self):
            return tuple([x for x in self._mu_T[:NUM_GENOTYPES]])
    
    property pi:
        def __get__(self):
            return tuple([x for x in self._pi[:NUM_JOINT_GENOTYPES]])
    
    cdef _normalise_pi(self):
        cdef int g
        cdef double norm_const
        
        norm_const = 0
        
        for g in range(NUM_JOINT_GENOTYPES):
            norm_const += self._pi[g]
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = self._pi[g] / norm_const

    cdef update(self, double * n, double * a_N, double * a_T, double * b_N, double * b_T):
        self._update_mu(self._mu_N, self._priors._mu_N, a_N, b_N)
        self._update_mu(self._mu_T, self._priors._mu_T, a_T, b_T)
        self._update_pi(n)

    cdef _update_mu(self, double * mu, double mu_prior[NUM_GENOTYPES][2], double * a, double * b):
        cdef int g
        cdef double alpha, beta, denom
    
        for g in range(NUM_GENOTYPES):
            alpha = a[g] + mu_prior[g][0] - 1
            beta = b[g] + mu_prior[g][1] - 1
            denom = alpha + beta

            mu[g] = alpha / denom
            
    cdef _update_pi(self, double * n):
        cdef int g
        cdef double denom
        cdef double pi[NUM_JOINT_GENOTYPES]
        
        denom = 0
        
        for g in range(NUM_JOINT_GENOTYPES):
            pi[g] = n[g] + self._priors._pi[g] - 1

        self._normalise_pi()
    
    cdef double _get_prior_log_likelihood(self):
        cdef double ll
        cdef double x[2]
        
        ll = 0
        
        for g in range(NUM_GENOTYPES):
            x[0] = self._mu_N[g]
            x[1] = 1 - self._mu_N[g]            
            ll += dirichlet_log_likelihood(x, self._priors._mu_N[g], NUM_BASES)
            
            x[0] = self._mu_T[g]
            x[1] = 1 - self._mu_T[g]
            ll += dirichlet_log_likelihood(x, self._priors._mu_T[g], NUM_BASES)
        
        ll += dirichlet_log_likelihood(self._pi, self._priors._pi, NUM_JOINT_GENOTYPES)
        
        return ll

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
