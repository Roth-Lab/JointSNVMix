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
        self._mu_N = params.mu_N
        self._mu_T = params.mu_T
        
        self._log_mix_weigths = [log(pi) for pi in self.params.pi]
        
        self._num_normal_genotypes = len(params.mu_N)
        
        self._num_tumour_genotypes = len(params.mu_T)
        
        self._num_joint_genotypes = self._num_normal_genotypes * self._num_tumour_genotypes

        self._init_ess()
    
    def _init_ess(self):
        '''
        Set the sufficient statistics to 0.
        '''
        self.a_N = [0] * self._num_normal_genotypes
        self.b_N = [0] * self._num_normal_genotypes
       
        self.a_T = [0] * self._num_tumour_genotypes
        self.b_T = [0] * self._num_tumour_genotypes
        
        self.n = [0] * self._num_joint_genotypes 
    
    def update(self, data_point, params):
        pass

class BinomialEss(JointSnvMixEss):
    def update(self, data_point, params):
        class_log_likelihood = [[0] * self._num_normal_genotypes] * self._num_tumour_genotypes
                
        for g_N in range(self._num_normal_genotypes):
            mu_N = self.params.mu_N[g_N]
            
            for g_T in range(self._num_tumour_genotypes):
                mu_T = self.params.mu_T[g_T]
                        
                class_log_likelihood[g_N][g_T] = binomial_log_likelihood(data_point.a_N, data_point.b_N, mu_N) + \
                                                 binomial_log_likelihood(data_point.a_T, data_point.b_T, mu_T)
            
        posterior = []
        
        for log_mix_weight, log_likelihood in zip(log_mix_weights, class_log_likelihood):
            posterior.append(log_mix_weight + log_likelihood)
        
        posterior = log_space_normalise(posterior)
        
        for resp in posterior:
            ess.n += resp
            
            ess.a_N += row[0] * resp
            ess.b_N += row[1] * resp
            ess.a_T += row[2] * resp
            ess.b_T += row[3] * resp
        
        return ess

#=======================================================================================================================
# CPT
#=======================================================================================================================
class Likelihood(object):
    def __init__(self, data, params):
        self._init_cpt_array(data, params)
        
        self._a_N = (< JointSnvMixOneData > data)._normal.counts[0]
        self._b_N = (< JointSnvMixOneData > data)._normal.counts[1]
        
        self._a_T = (< JointSnvMixOneData > data)._tumour.counts[0]
        self._b_T = (< JointSnvMixOneData > data)._tumour.counts[1]

    cdef double * get_resp(self):
        cdef int g
        cdef double * resp
    
        resp = < double *> malloc(NUM_JOINT_GENOTYPES * sizeof(double))
        
        for g in range(NUM_JOINT_GENOTYPES):
            resp[g] = self._cpt_array[g]
        
        log_space_normalise_row(resp, NUM_JOINT_GENOTYPES)
        
        for g in range(NUM_JOINT_GENOTYPES):
            resp[g] = exp(resp[g])
        
        return resp
    
    cdef double * get_expected_counts_a_N(self):
        cdef int g
        cdef double * marginal_resp, * a
        
        marginal_resp = self._get_normal_marginal_resp()
        
        a = self._get_expected_counts(self._a_N, marginal_resp)
        
        free(marginal_resp)
        
        return a

    cdef double * get_expected_counts_a_T(self):
        cdef int g
        cdef double * marginal_resp, * a
        
        marginal_resp = self._get_tumour_marginal_resp()
        
        a = self._get_expected_counts(self._a_T, marginal_resp)
        
        free(marginal_resp)
        
        return a
    
    cdef double * get_expected_counts_b_N(self):
        cdef int g
        cdef double * marginal_resp, * b
        
        marginal_resp = self._get_normal_marginal_resp()
        
        b = self._get_expected_counts(self._b_N, marginal_resp)
        
        free(marginal_resp)
        
        return b

    cdef double * get_expected_counts_b_T(self):
        cdef int g
        cdef double * marginal_resp, * b
        
        marginal_resp = self._get_tumour_marginal_resp()
        
        b = self._get_expected_counts(self._b_T, marginal_resp)
        
        free(marginal_resp)
        
        return b

    cdef double get_log_sum(self):
        cdef int g
        cdef double log_marginal
        
        log_marginal = 0
        
        log_marginal = log_sum_exp(& self._cpt_array[0], NUM_JOINT_GENOTYPES)
        
        return log_marginal

    cdef double * _get_normal_marginal_resp(self):
        cdef int g_N, g_T, g_J
        cdef double * resp, * marginal_resp
    
        resp = self.get_resp()
        
        marginal_resp = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g_N in range(NUM_GENOTYPES):
            marginal_resp[g_N] = 0
            for g_T in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                marginal_resp[g_N] += resp[g_J]
        
        free(resp)
        
        return marginal_resp

    cdef double * _get_tumour_marginal_resp(self):
        cdef int g_N, g_T, g_J
        cdef double * resp, * marginal_resp
    
        resp = self.get_resp()
        
        marginal_resp = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g_T in range(NUM_GENOTYPES):
            marginal_resp[g_T] = 0
            for g_N in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                marginal_resp[g_T] += resp[g_J]
        
        free(resp)
        
        return marginal_resp

    cdef double * _get_expected_counts(self, int counts, double * marginal_resp):
        cdef int g
        cdef double * expected_counts
        
        expected_counts = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            expected_counts[g] = counts * marginal_resp[g]
        
        return expected_counts
    
    cdef _init_cpt_array(self, JointSnvMixOneData data, JointSnvMixParameters params):
        cdef int a_N, b_N, a_T, b_T, g_N, g_T, g_J
        cdef double mu_N, mu_T, log_pi

        a_N = data._normal.counts[0]
        b_N = data._normal.counts[1]
        a_T = data._tumour.counts[0]
        b_T = data._tumour.counts[1]   
        
        for g_N in range(NUM_GENOTYPES):      
            mu_N = params._mu_N[g_N]
            
            for g_T in range(NUM_GENOTYPES):         
                mu_T = params._mu_T[g_T]
            
                g_J = NUM_GENOTYPES * g_N + g_T
            
                log_pi = log(params._pi[g_J])
            
                self._cpt_array[g_J] = log_pi + \
                                       self._binomial_log_likelihood(a_N, b_N, mu_N) + \
                                       self._binomial_log_likelihood(a_T, b_T, mu_T)
        
    cdef double _binomial_log_likelihood(self, int a, int b, double mu):
        return a * log(mu) + b * log(1 - mu)
