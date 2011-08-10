'''
Created on 2011-08-04

@author: Andrew Roth
'''
DEF FLOAT_INFN = float('-inf')

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9
DEF NUM_BASES = 2
DEF EPS = 1e-100

#=======================================================================================================================
# Subsamplers
#=======================================================================================================================
cdef class PairedDataSubSampler(object):   
    def __init__(self, int skip_size, int min_normal_depth, int min_tumour_depth):
        self._skip_size = skip_size
        self._min_normal_depth = min_normal_depth
        self._min_tumour_depth = min_tumour_depth

    def subsample(self, Counter counter, refs=None):
        cdef int i, ref_sample_size
        cdef RefIterator ref_iter       
        cdef PairedSampleBinomialCounterRow row
        cdef list sample
        
        sample = []
        
        print '''Randomly sub-sampling every {0}th position with normal depth {1} and tumour depth {2} the data set.'''.format(self._skip_size, self._min_normal_depth, self._min_tumour_depth)
        
        if refs == None:
            refs = counter.refs
        
        for ref in refs:
            print "Subsampling ref {0}.".format(ref)

            ref_iter = counter.iter_ref(ref)
            
            i = 0
            ref_sample_size = 0
            
            try:
                while True:
                    ref_iter.cnext()
                    
                    row = ref_iter._current_row
                    
                    if row._normal_depth < self._min_normal_depth or row._tumour_depth < self._min_tumour_depth:
                        continue
                    
                    if i % self._skip_size == 0:               
                        self._add_row_to_sample(sample, row)
                        ref_sample_size += 1
                        
                    i += 1
                        
            except StopIteration:
                pass
            
            print "Sub-sampled {0} positions from ref {1}".format(ref_sample_size, ref)
        
        print "Total sub-sample size is {0}".format(len(sample))
        
        return sample
    
    cdef _add_row_to_sample(self, list sample, PairedSampleBinomialCounterRow row):
        pass
 
cdef class JointSnvMixOneSubsampler(PairedDataSubSampler):
    cdef _add_row_to_sample(self, list sample, PairedSampleBinomialCounterRow row):
        cdef JointSnvMixOneData data
    
        data = makeJointSnvMixOneData(row)
        
        sample.append(data)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixTwoSubsampler(PairedDataSubSampler):
    cdef _add_row_to_sample(self, list sample, PairedSampleBinomialCounterRow row):
        cdef JointSnvMixTwoData data
    
        data = makeJointSnvMixTwoData(row)
        
        sample.append(data)

#=======================================================================================================================
# Data
#=======================================================================================================================
cdef class JointSnvMixData(object):
    pass
 
cdef class JointSnvMixOneData(JointSnvMixData):
    pass

cdef JointSnvMixOneData makeJointSnvMixOneData(JointBinaryCounterRow row):
    cdef JointSnvMixOneData data = JointSnvMixOneData.__new__(JointSnvMixOneData)
    
    data._normal = makeSnvMixOneData(row._normal_counts)
    data._tumour = makeSnvMixOneData(row._tumour_counts)
    
    return data

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixTwoData(JointSnvMixData):
    pass

cdef JointSnvMixTwoData makeJointSnvMixTwoData(JointBinaryQualityCounterRow row):
    cdef JointSnvMixTwoData data = JointSnvMixTwoData.__new__(JointSnvMixTwoData)
    
    data._normal = makeSnvMixTwoData(row._normal_data)
    data._tumour = makeSnvMixTwoData(row._tumour_data)
    
    return data

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
        
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        
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
        s = "mu_N_alpha : {0}, {1}, {2}\n".format(self.mu_N[0][0],
                                                  self.mu_N[1][0],
                                                  self.mu_N[2][0])
        
        s += "mu_N_beta : {0}, {1}, {2}\n".format(self.mu_N[0][1],
                                                  self.mu_N[1][1],
                                                  self.mu_N[2][1])
        
        s = "mu_T_alpha : {0}, {1}, {2}\n".format(self.mu_T[0][0],
                                                  self.mu_T[1][0],
                                                  self.mu_T[2][0])
        
        s += "mu_T_beta : {0}, {1}, {2}\n".format(self.mu_T[0][1],
                                                  self.mu_T[1][1],
                                                  self.mu_T[2][1])
        
        s += "pi : "
        s += "\t".join([str(x) for x in self.pi[:NUM_JOINT_GENOTYPES]])
        s += "\n"
        
        return s


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
        
        norm_const = sum(pi)
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = pi[g] / norm_const
            
        
    def __str__(self):
        s = "mu_N : "
        s += "\t".join([str(x) for x in self._mu_N[:NUM_GENOTYPES]])
        s += "\n"
        
        s += "mu_T : "
        s += "\t".join([str(x) for x in self._mu_T[:NUM_GENOTYPES]])
        s += "\n"

        s += "pi : "
        s += "\t".join([str(x) for x in self._pi[:NUM_JOINT_GENOTYPES]])
        s += "\n"
        
        return s
    
    property mu_N:
        def __get__(self):
            return tuple([x for x in self._mu_N[:NUM_GENOTYPES]])
    
    property mu_T:
        def __get__(self):
            return tuple([x for x in self._mu_T[:NUM_GENOTYPES]])
    
    property pi:
        def __get__(self):
            return tuple([x for x in self._pi[:NUM_JOINT_GENOTYPES]])

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
            denom += pi[g]
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = pi[g] / denom
    
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
# Models
#=======================================================================================================================       
cdef class JointSnvMixModel(object):
    def __init__(self, JointSnvMixParameters params):
        self._params = params
    
    def fit(self, list data, convergence_threshold=1e-6, max_iters=100):
        trainer = JointSnvMixModelTrainer(self, convergence_threshold, max_iters)
        
        trainer.train(data)
        
    property params:
        def __get__(self):
            return self._params

    cdef double _get_lower_bound(self, list data):
        cdef double lb
        cdef JointSnvMixData pos_data
        
        lb = 0
        
        for pos_data in data:
            lb += self._get_log_likelihood(pos_data)
        
        lb += self._params._get_prior_log_likelihood()
        
        return lb
        
    cdef JointSnvMixCpt _get_complete_log_likelihood(self, JointSnvMixData data):
        pass
    
    cdef double _get_log_likelihood(self, JointSnvMixData data):
        cdef JointSnvMixCpt cpt
        cdef double likelihood
        
        cpt = self._get_complete_log_likelihood(data)
        
        likelihood = cpt.marginalise()
        
        if likelihood == 0:
            likelihood = EPS
        
        return log(likelihood)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixOneModel(JointSnvMixModel):
    cdef JointSnvMixCpt _get_complete_log_likelihood(self, JointSnvMixData data):
        return JointSnvMixOneCpt(data, self._params)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixTwoModel(JointSnvMixModel):
    cdef JointSnvMixCpt _get_complete_log_likelihood(self, JointSnvMixData data):
        return JointSnvMixTwoCpt(data, self._params)
    
#=======================================================================================================================
# Trainer
#=======================================================================================================================
cdef class JointSnvMixModelTrainer(object):    
    def __init__(self, model, convergence_threshold, max_iters):
        self._model = model
        self._convergence_threshold = convergence_threshold
        self._max_iters = max_iters
        
        self._converged = 0
        self._iters = 0
        self._lower_bounds = [FLOAT_INFN]

    cdef train(self, list data):
        cdef JointSnvMixEss ess
        
        while not self._converged:
            ess = self._do_e_step(data)
            self._do_m_step(ess)
            
            self._check_convergence(data)
            
            print self._iters, self._lower_bounds[-1]
            print self._model.params

    cdef JointSnvMixEss _do_e_step(self, list data):                              
        cdef JointSnvMixData pos_data
        cdef JointSnvMixEss ess
        cdef JointSnvMixCpt cpt
        
        ess = JointSnvMixEss()
    
        for pos_data in data:
            cpt = self._model._get_complete_log_likelihood(pos_data)          
        
            ess.update(cpt)
        
        return ess
    
    cdef void _do_m_step(self, JointSnvMixEss ess):       
        self._model._params.update(ess._n, ess._a_N, ess._a_T, ess._b_N, ess._b_T)

    cdef _check_convergence(self, list data):
        cdef double rel_change, lb, ll, prev_ll
        
        lb = self._model._get_lower_bound(data)        
        self._lower_bounds.append(lb)
        
        ll = self._lower_bounds[-1]
        prev_ll = self._lower_bounds[-2]
        
        rel_change = (ll - prev_ll) / abs(prev_ll)
    
        if rel_change < 0:
            print "Lower bound decreased exiting."
            self._converged = 1
        elif rel_change < self._convergence_threshold:
            print "Converged"
            self._converged = 1
        elif self._iters >= self._max_iters:
            print "Maximum number of iters exceeded exiting."
            self._converged = 1
        else:
            self._converged = 0
        
        self._iters += 1

#=======================================================================================================================
# ESS
#=======================================================================================================================
cdef class JointSnvMixEss(object):
    def __init__(self):
        self.reset()
    
    cdef void reset(self):
        cdef int g
        
        for g in range(NUM_GENOTYPES):
            self._a_N[g] = 0
            self._a_T[g] = 0
            self._b_N[g] = 0
            self._b_T[g] = 0
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._n[g] = 0
    
    cdef update(self, JointSnvMixCpt cpt):
        cdef double * a_N, * a_T, * b_N, * b_T, * resp
        
        resp = cpt.get_resp()
        a_N = cpt.get_expected_counts_a_N()
        a_T = cpt.get_expected_counts_a_T()
        b_N = cpt.get_expected_counts_b_N()
        b_T = cpt.get_expected_counts_b_T()
    
        for g in range(NUM_GENOTYPES):            
            self._a_N[g] += a_N[g]
            self._a_T[g] += a_T[g]
            
            self._b_N[g] += b_N[g]
            self._b_T[g] += b_T[g]
            
            print a_N[g], a_T[g], b_N[g], b_T[g]
        
        print
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._n[g] += resp[g]
            
            print resp[g],
            
        print
        
        free(resp)
        free(a_N)
        free(a_T)
        free(b_N)
        free(b_T)

#=======================================================================================================================
# CPT
#=======================================================================================================================
cdef class JointSnvMixCpt(object):
    cdef double * get_resp(self):
        pass
    
    cdef double * get_expected_counts_a_N(self):
        pass

    cdef double * get_expected_counts_a_T(self):
        pass
    
    cdef double * get_expected_counts_b_N(self):
        pass

    cdef double * get_expected_counts_b_T(self):
        pass
    
    cdef double marginalise(self):
        pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixOneCpt(JointSnvMixCpt):
    def __init__(self, JointSnvMixData data, JointSnvMixParameters params):
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

    cdef double marginalise(self):
        cdef int g
        cdef double marginal
        
        marginal = 0
        
        for g in range(NUM_JOINT_GENOTYPES):
            marginal += exp(self._cpt_array[g])
        
        return marginal

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
        
        for g_N in range(NUM_GENOTYPES):
            a_N = data._normal.counts[0]
            b_N = data._normal.counts[1]            
            mu_N = params.mu_N[g_N]
            
            for g_T in range(NUM_GENOTYPES):
                a_T = data._tumour.counts[0]
                b_T = data._tumour.counts[1]            
                mu_T = params.mu_T[g_T]
            
                g_J = NUM_GENOTYPES * g_N + g_T
            
                log_pi = log(params.pi[g_J])
            
                self._cpt_array[g_J] = log_pi + \
                                       self._binomial_log_likelihood(a_N, b_N, mu_N) + \
                                       self._binomial_log_likelihood(a_T, b_T, mu_T)
        
    cdef double _binomial_log_likelihood(self, int a, int b, double mu):
        return a * log(mu) + b * log(1 - mu)
        
#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixTwoCpt(JointSnvMixCpt):
    def __init__(self, JointSnvMixTwoData data, JointSnvMixParameters params):
        self._init_cpt(data, params)

    cdef double * get_resp(self):
        cdef int g
        cdef double * resp
        
        resp = < double *> malloc(NUM_JOINT_GENOTYPES * sizeof(double))
        
        for g in range(NUM_JOINT_GENOTYPES):     
            resp[g] = self._resp[g]
        
        return resp
    
    cdef double * get_expected_counts_a_N(self):
        cdef int g
        cdef double * counts
        
        counts = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):     
            counts[g] = self._normal_counts_a[g]
        
        return counts
    
    cdef double * get_expected_counts_a_T(self):
        cdef int g
        cdef double * counts
        
        counts = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):     
            counts[g] = self._tumour_counts_a[g]
        
        return counts
    cdef double * get_expected_counts_b_N(self):
        cdef int g
        cdef double * counts
        
        counts = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):     
            counts[g] = self._normal_counts_b[g]
        
        return counts
    
    cdef double * get_expected_counts_b_T(self):
        cdef int g
        cdef double * counts
        
        counts = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):     
            counts[g] = self._tumour_counts_b[g]
        
        return counts
    
    cdef double marginalise(self):
        return self._marginal
    
    cdef void _init_cpt(self, JointSnvMixTwoData data, JointSnvMixParameters params):
        cdef double **** normal_cpt, **** tumour_cpt
        cdef double ** normal_read_marginals, ** tumour_read_marginals
        cdef double * normal_class_marginals, * tumour_class_marginals, * joint_class_marginals
        
        normal_cpt = self._get_cpt_array(data._normal, params._mu_N)
        tumour_cpt = self._get_cpt_array(data._tumour, params._mu_T)
        
        normal_read_marginals = self._get_read_marginals(normal_cpt, data._normal.depth)
        tumour_read_marginals = self._get_read_marginals(tumour_cpt, data._tumour.depth)
        
        normal_class_marginals = self._get_class_marginals(normal_read_marginals, data._normal.depth)
        tumour_class_marginals = self._get_class_marginals(tumour_read_marginals, data._tumour.depth)
        
        joint_class_marginals = self._get_joint_class_marginals(normal_class_marginals,
                                                                tumour_class_marginals,
                                                                params._pi)
                
        self._init_marginal(normal_class_marginals, tumour_class_marginals, params._pi)
        self._init_resp(normal_class_marginals, tumour_class_marginals, params._pi)
        self._init_normal_expected_counts(normal_cpt, normal_read_marginals, normal_class_marginals, data._normal.depth)
        self._init_tumour_expected_counts(tumour_cpt, tumour_read_marginals, tumour_class_marginals, data._tumour.depth)
        
        self._free_cpt_array(normal_cpt, data._normal.depth)
        self._free_cpt_array(tumour_cpt, data._tumour.depth)
        
        self._free_read_marginals(normal_read_marginals)
        self._free_read_marginals(tumour_read_marginals)
        
        free(normal_class_marginals)
        free(tumour_class_marginals)
        free(joint_class_marginals)

    cdef double **** _get_cpt_array(self, SnvMixTwoData data, double * mu):
        cdef int g, d, a, z
        cdef double m, q, r 

        cpt_array = self._make_cpt_array(data.depth)

        for g in range(NUM_GENOTYPES):
            m = mu[g]
            
            for d in range(data.depth):   
                q = data.q[d]
                r = data.r[d]
                
                for a in range(2):
                    for z in range(2):                        
                        cpt_array[g][d][a][z] = self._get_read_complete_likelihood(a, z, q, r, m)
        
        return cpt_array
    
    cdef double ** _get_read_marginals(self, double **** cpt_array, int depth):
        cdef int  g, d, a, z
        cdef double ** read_marginals
        cdef double temp_marginal[4]
        
        read_marginals = self._make_read_marginals(depth)

        for g in range(NUM_GENOTYPES):
            for d in range(depth):
                for a in range(2):
                    for z in range(2):                        
                        temp_marginal[a * 2 + z] = cpt_array[g][d][a][z]
                
                read_marginals[g][d] = log_sum_exp(temp_marginal, 4)

        return read_marginals
    
    cdef double * _get_class_marginals(self, double ** read_marginals, int depth):
        cdef int  g, d, a, z
        cdef double * class_marginals
        
        class_marginals = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            class_marginals[g] = 0
            
            for d in range(depth):
                class_marginals[g] += read_marginals[g][d]
        
        return class_marginals

    cdef double * _get_joint_class_marginals(self, double * normal_marginals, double * tumour_marginals, double * pi):
        cdef int g_N, g_T, g_J
        cdef double * joint_marginals
        
        joint_marginals = < double *> malloc(NUM_JOINT_GENOTYPES * sizeof(double))
        
        for g_N in range(NUM_GENOTYPES):
            for g_T in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                joint_marginals[g_J] = log(pi[g_J]) + normal_marginals[g_N] + tumour_marginals[g_T]
    
        return joint_marginals
    
    cdef void _init_marginal(self, double * normal_marginals, double * tumour_marginals, double * pi):
        cdef int g_N, g_T, g_J
        cdef double log_norm_const
        cdef double ll[NUM_JOINT_GENOTYPES]
        
        for g_N in range(NUM_GENOTYPES):
            for g_T in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                ll[g_J] = log(pi[g_J]) + normal_marginals[g_N] + tumour_marginals[g_T]
        
        norm_const = log_sum_exp(ll, NUM_JOINT_GENOTYPES)
        norm_const = exp(norm_const)
        
        if norm_const == 0:
            norm_const = EPS
        
        self._marginal = norm_const
    
    cdef void _init_resp(self, double * normal_marginals, double * tumour_marginals, double * pi):
        cdef int g_N, g_T, g_J
        
        for g_N in range(NUM_GENOTYPES):
            for g_T in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                self._resp[g_J] = log(pi[g_J]) + normal_marginals[g_N] + tumour_marginals[g_T]
                
        log_space_normalise_row(self._resp, NUM_JOINT_GENOTYPES)
        
        for g_J in range(NUM_JOINT_GENOTYPES):
            self._resp[g_J] = exp(self._resp[g_J])

    cdef void _init_normal_expected_counts(self,
                                           double **** cpt_array,
                                           double ** read_marginals,
                                           double * joint_marginals,
                                           int depth):
        cdef int g_N, g_T, d
        cdef double norm_const[NUM_GENOTYPES]
        cdef double temp[NUM_GENOTYPES]
        cdef double x[2]
        
        for g_N in range(NUM_GENOTYPES):
            norm_const[g_N] = 0
            for g_T in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                temp[g_T] = joint_marginals[g_J]
                        
            norm_const[g_N] = log_sum_exp(temp, NUM_GENOTYPES)
            
            norm_const[g_N] = norm_const[g_N] - log(self._marginal) 

        for g_N in range(NUM_GENOTYPES):
            self._normal_counts_a[g_N] = 0
            self._normal_counts_b[g_N] = 0
            
            for d in range(depth):   
                read_prob = norm_const[g_N] - read_marginals[g_N][d]
                
                x[0] = cpt_array[g_N][d][1][0]
                x[1] = cpt_array[g_N][d][1][1]
                
                self._normal_counts_a[g_N] += exp(read_prob + log_sum_exp(x, 2))
                
                x[0] = cpt_array[g_N][d][0][0]
                x[1] = cpt_array[g_N][d][0][1]
                self._normal_counts_b[g_N] += exp(read_prob + log_sum_exp(x, 2))

    cdef void _init_tumour_expected_counts(self,
                                           double **** cpt_array,
                                           double ** read_marginals,
                                           double * joint_marginals,
                                           int depth):
        cdef int g_N, g_T, d
        cdef double norm_const[NUM_GENOTYPES]
        cdef double temp[NUM_GENOTYPES]
        cdef double x[2]
        
        for g_T in range(NUM_GENOTYPES):
            norm_const[g_N] = 0
            for g_N in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                temp[g_N] = joint_marginals[g_J]
                        
            norm_const[g_T] = log_sum_exp(temp, NUM_GENOTYPES)
            
            norm_const[g_T] = norm_const[g_T] - log(self._marginal) 

        for g_T in range(NUM_GENOTYPES):
            self._tumour_counts_a[g_T] = 0
            self._tumour_counts_b[g_T] = 0
            
            for d in range(depth):
                read_prob = norm_const[g_T] - read_marginals[g_T][d]
                
                x[0] = cpt_array[g_T][d][1][0]
                x[1] = cpt_array[g_T][d][1][1]
                
                self._tumour_counts_a[g_T] += exp(read_prob + log_sum_exp(x, 2))
                
                x[0] = cpt_array[g_T][d][0][0]
                x[1] = cpt_array[g_T][d][0][1]
                self._tumour_counts_b[g_T] += exp(read_prob + log_sum_exp(x, 2))

    cdef double _get_read_complete_likelihood(self, int a, int z, double q, double r, double mu):
        if a == 0 and z == 0:
            return log(0.5) + log(1 - r) + log(1 - mu)
        elif a == 0 and z == 1:
            return log(1 - q) + log(r) + log(1 - mu)
        elif a == 1 and z == 0:
            return log(0.5) + log(1 - r) + log(mu)
        else:
            return log(q) + log(r) + log(mu)
#    cdef double _get_read_complete_likelihood(self, int a, int z, double q, double r, double mu):
#        if a == 0 and z == 0:
#            return 0.5 * (1 - r) * (1 - mu)
#        elif a == 0 and z == 1:
#            return (1 - q) * r * (1 - mu)
#        elif a == 1 and z == 0:
#            return 0.5 * (1 - r) * mu
#        else:
#            return q * r * mu        

    cdef double **** _make_cpt_array(self, int depth):
        cdef int g, d, a, z
        
        cpt_array = < double ****> malloc(NUM_GENOTYPES * sizeof(double *))
        
        for g in range(NUM_GENOTYPES):
            cpt_array[g] = < double ***> malloc(depth * sizeof(double *))
            
            for d in range(depth):
                cpt_array[g][d] = < double **> malloc(2 * sizeof(double *))
                
                for a in range(2):
                    cpt_array[g][d][a] = < double *> malloc(2 * sizeof(double))
                    
                    for z in range(2):
                        cpt_array[g][d][a][z] = 0
        
        return cpt_array

    cdef void _free_cpt_array(self, double **** cpt_array, int depth):
        cdef int d, g, a
        
        for g in range(NUM_GENOTYPES):
            for d in range(depth):
                for a in range(2):
                    free(cpt_array[g][d][a])
                free(cpt_array[g][d])        
            free(cpt_array[g])
        free(cpt_array)
    
    cdef double ** _make_read_marginals(self, int depth):
        cdef double ** read_marigals
        
        read_marginals = < double **> malloc(NUM_GENOTYPES * sizeof(double))

        for g in range(NUM_GENOTYPES):
            read_marginals[g] = < double *> malloc(depth * sizeof(double))
        
        return read_marginals
        

    cdef void _free_read_marginals(self, double ** read_marginals):
        for g in range(NUM_GENOTYPES):
            free(read_marginals[g])
        
        free(read_marginals)
