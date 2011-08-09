'''
Created on 2011-08-04

@author: Andrew Roth
'''
DEF FLOAT_INFN = float('-inf')

DEF NUM_GENOTYPES = 3
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
        cdef dict sample
        
        sample = {'normal' : [], 'tumour' : []}
        
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
        
        print "Total sub-sample size is {0}".format(len(sample['normal']))
        
        return sample
    
    cdef _add_row_to_sample(self, dict sample, PairedSampleBinomialCounterRow row):
        pass
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneSubsampler(PairedDataSubSampler):
    cdef _add_row_to_sample(self, dict sample, PairedSampleBinomialCounterRow row):
        cdef SnvMixOneData normal_data, tumour_data
    
        normal_data = makeSnvMixOneData((< JointBinaryCounterRow > row)._normal_counts)
        tumour_data = makeSnvMixOneData((< JointBinaryCounterRow > row)._tumour_counts)
        
        sample['normal'].append(normal_data)
        sample['tumour'].append(tumour_data)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoSubsampler(PairedDataSubSampler):
    cdef _add_row_to_sample(self, dict sample, PairedSampleBinomialCounterRow row):
        cdef SnvMixTwoData normal_data, tumour_data
    
        normal_data = makeSnvMixTwoData((< JointBinaryQualityCounterRow > row)._normal_data)
        tumour_data = makeSnvMixTwoData((< JointBinaryQualityCounterRow > row)._tumour_data)
        
        sample['normal'].append(normal_data)
        sample['tumour'].append(tumour_data)


#=======================================================================================================================
# Data
#=======================================================================================================================
cdef class SnvMixData(object):
    pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneData(SnvMixData):
    pass

cdef SnvMixOneData makeSnvMixOneData(binary_counts_struct counts):
    cdef SnvMixOneData data = SnvMixOneData.__new__(SnvMixOneData)
    
    data.counts[0] = counts.A
    data.counts[1] = counts.B
    
    return data

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoData(SnvMixData):
    def __dealloc__(self):
        free(self.labels)
        free(self.q)
        free(self.r)

cdef SnvMixTwoData makeSnvMixTwoData(base_map_qualities_struct data_struct):
    cdef read_index, i, l
    cdef double temp_q
    cdef SnvMixTwoData data = SnvMixTwoData.__new__(SnvMixTwoData)
    
    l = data_struct.depth.A + data_struct.depth.B
    
    data.depth = l
    
    data.labels = < int *> malloc(l * sizeof(int))
    data.q = < double *> malloc(l * sizeof(double))
    data.r = < double *> malloc(l * sizeof(double))
    
    i = 0
    
    for read_index in range(data_struct.depth.A):
        data.labels[i] = 1
        data.q[i] = get_phred_qual_to_prob(data_struct.base_quals.A[read_index])
        data.r[i] = get_phred_qual_to_prob(data_struct.map_quals.A[read_index])
        
        i += 1
    
    for read_index in range(data_struct.depth.B):
        data.labels[i] = 0
        
        temp_q = get_phred_qual_to_prob(data_struct.base_quals.B[read_index])
        data.q[i] = (1 - temp_q) / 3
                
        data.r[i] = get_phred_qual_to_prob(data_struct.map_quals.B[read_index])
        
        i += 1
    
    return data

cdef double get_phred_qual_to_prob(int qual):
    cdef double base, exp, prob
    
    exp = -1 * (< double > qual) / 10
    base = 10
    
    prob = 1 - pow(base, exp)

    return prob

#=======================================================================================================================
# Priors
#=======================================================================================================================
cdef class SnvMixPriors(object):
    def __init__(self, **kwargs):
        default_mu = (
                      (100, 2),
                      (50, 50),
                      (2, 100)
                      )
        
        default_pi = (1000, 100, 100)
        
        mu = kwargs.get('mu', default_mu)
        pi = kwargs.get('pi', default_pi)
        
        for g in range(NUM_GENOTYPES):
            self.mu[g][0] = mu[g][0]
            self.mu[g][1] = mu[g][1]
            self.pi[g] = pi[g]

    def __str__(self):
        s = "mu_alpha : {0}, {1}, {2}\n".format(self.mu[0][0],
                                                self.mu[1][0],
                                                self.mu[2][0])
        
        s += "mu_beta : {0}, {1}, {2}\n".format(self.mu[0][1],
                                                self.mu[1][1],
                                                self.mu[2][1])
        
        s += "pi : {0}, {1}, {2}\n".format(self.pi[0],
                                           self.pi[1],
                                           self.pi[2])
        
        return s

#---------------------------------------------------------------------------------------------------------------------- 
cdef class PairedSnvMixPriors(object):
    def __init__(self, **kwargs):
        default_mu = (
                      (100, 2),
                      (50, 50),
                      (2, 100)
                      )
        
        default_pi = (1000, 100, 100)
        
        mu_N = kwargs.get('mu_N', default_mu)
        mu_T = kwargs.get('mu_T', default_mu)

        pi_N = kwargs.get('pi_N', default_pi)
        pi_T = kwargs.get('pi_T', default_pi)
        
        self._normal_priors = SnvMixPriors(mu=mu_N, pi=pi_N)
        self._tumour_priors = SnvMixPriors(mu=mu_T, pi=pi_T)
    
    property normal:
        def __get__(self):
            return self._normal_priors
    
    property tumour:
        def __get__(self):
            return self._tumour_priors


#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class SnvMixParameters(object):
    def __init__(self, **kwargs):        
        self.priors = kwargs.get('priors', SnvMixPriors())
        
        mu = kwargs.get('mu', (0.99, 0.5, 0.01))
        pi = kwargs.get('pi', (0.99, 0.009, 0.001))
        
        for g in range(NUM_GENOTYPES):
            self.mu[g] = mu[g]
            self.pi[g] = pi[g]
        
    def __str__(self):
        s = "mu : {0}, {1}, {2}\n".format(self.mu[0],
                                          self.mu[1],
                                          self.mu[2])
        
        s += "pi : {0}, {1}, {2}\n".format(self.pi[0],
                                           self.pi[1],
                                           self.pi[2])
        
        return s

    cdef update(self, double * n, double * a, double * b):
        self._update_mu(a, b)
        self._update_pi(n)

    cdef _update_mu(self, double * a, double * b):
        cdef int g
        cdef double alpha, beta, denom
    
        for g in range(NUM_GENOTYPES):
            alpha = a[g] + self.priors.mu[g][0] - 1
            beta = b[g] + self.priors.mu[g][1] - 1
            denom = alpha + beta

            self.mu[g] = alpha / denom
            
    cdef _update_pi(self, double * n):
        cdef int g
        cdef double denom
        cdef double pi[NUM_GENOTYPES]
        
        denom = 0
        
        for g in range(NUM_GENOTYPES):
            pi[g] = n[g] + self.priors.pi[g] - 1
            denom += pi[g]
        
        for g in range(NUM_GENOTYPES):
            self.pi[g] = pi[g] / denom
    
    cdef double _get_prior_log_likelihood(self):
        cdef double ll
        cdef double x[2]
        
        ll = 0
        
        for g in range(NUM_GENOTYPES):
            x[0] = self.mu[g]
            x[1] = 1 - self.mu[g]
            
            ll += dirichlet_log_likelihood(x, self.priors.mu[g], NUM_BASES)
        
        ll += dirichlet_log_likelihood(self.pi, self.priors.pi, NUM_GENOTYPES)
        
        return ll
            
#---------------------------------------------------------------------------------------------------------------------- 
cdef class PairedSnvMixParameters(object):
    def __init__(self, **kwargs):
        priors = kwargs.get('priors', PairedSnvMixPriors())
         
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        self._normal_params = SnvMixParameters(priors=priors._normal_priors, mu=mu_N, pi=pi_N)
        self._tumour_params = SnvMixParameters(priors=priors._tumour_priors, mu=mu_T, pi=pi_T)
    
    property normal:
        def __get__(self):
            return self._normal_params
    
    property tumour:
        def __get__(self):
            return self._tumour_params

#=======================================================================================================================
# Models
#=======================================================================================================================       
cdef class SnvMixModel(object):
    def __init__(self, SnvMixParameters params):
        self.params = params
    
    def fit(self, list data, convergence_threshold=1e-6, max_iters=100):
        trainer = SnvMixModelTrainer(self, convergence_threshold, max_iters)
        
        trainer.train(data)
    
    def predict(self, SnvMixData data):
        pass

    cdef double _get_lower_bound(self, list data):
        cdef double lb
        cdef SnvMixData pos_data
        
        lb = 0
        
        for pos_data in data:
            lb += self._get_log_likelihood(pos_data)
        
        lb += self.params._get_prior_log_likelihood()
        
        return lb
        
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data):
        pass
    
    cdef double _get_log_likelihood(self, SnvMixData data):
        cdef SnvMixCpt cpt
        cdef double likelihood
        
        cpt = self._get_complete_log_likelihood(data)
        
        likelihood = cpt.marginalise()
        
        if likelihood == 0:
            likelihood = EPS
        
        return log(likelihood)

cdef class PairedSnvMixModel(object):
    def __init__(self, SnvMixModel normal_model, SnvMixModel tumour_model):
        self.normal_model = normal_model
        self.tumour_model = tumour_model
        
    def fit(self, list normal_data, list tumour_data):
        self.normal_model.fit(normal_data)
        self.tumour_model.fit(tumour_data)
    
    def predict(self, SnvMixData normal_data, SnvMixData tumour_data):
        normal_resp = self.normal_model.predict(normal_data)
        tumour_resp = self.tumour_model.predict(tumour_data)
        
        joint_resp = self._get_joint_resp(normal_resp, tumour_resp)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneModel(SnvMixModel):
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data):
        return SnvMixOneCpt(data, self.params)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoModel(SnvMixModel):
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data):
        return SnvMixTwoCpt(data, self.params)
    
#=======================================================================================================================
# Trainer
#=======================================================================================================================
cdef class SnvMixModelTrainer(object):    
    def __init__(self, model, convergence_threshold, max_iters):
        self._model = model
        self._convergence_threshold = convergence_threshold
        self._max_iters = max_iters
        
        self._converged = 0
        self._iters = 0
        self._lower_bounds = [FLOAT_INFN]

    cdef train(self, list data):
        cdef SnvMixEss ess
        
        while not self._converged:
            ess = self._do_e_step(data)
            self._do_m_step(ess)
            
            self._check_convergence(data)
            
            print self._iters, self._lower_bounds[-1]
            print self._model.params

    cdef SnvMixEss _do_e_step(self, list data):                              
        cdef SnvMixData pos_data
        cdef SnvMixEss ess
        cdef SnvMixCpt cpt
        
        ess = SnvMixEss()
    
        for pos_data in data:
            cpt = self._model._get_complete_log_likelihood(pos_data)                   
        
            ess.update(cpt)
        
        return ess
    
    cdef void _do_m_step(self, SnvMixEss ess):       
        self._model.params.update(ess.n, ess.a, ess.b)

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
cdef class SnvMixEss(object):
    def __init__(self):
        self.reset()
    
    cdef void reset(self):
        cdef int g
        
        for g in range(NUM_GENOTYPES):
            self.a[g] = 0
            self.b[g] = 0
            self.n[g] = 0
    
    cdef update(self, SnvMixCpt cpt):
        cdef double * a, *b, * resp
        
        resp = cpt.get_resp()
        a = cpt.get_expected_counts_a()
        b = cpt.get_expected_counts_b()
    
        for g in range(NUM_GENOTYPES):
            self.n[g] += resp[g]
            self.a[g] += a[g]
            self.b[g] += b[g]
        
        free(resp)
        free(a)
        free(b)

#=======================================================================================================================
# CPT
#=======================================================================================================================
cdef class SnvMixCpt(object):
    cdef double * get_resp(self):
        pass
    
    cdef double * get_expected_counts_a(self):
        pass
    
    cdef double * get_expected_counts_b(self):
        pass
    
    cdef double marginalise(self):
        pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneCpt(SnvMixCpt):
    def __init__(self, SnvMixData data, SnvMixParameters params):
        self._init_cpt_array(data, params)
        
        self._a = (< SnvMixOneData > data).counts[0]
        self._b = (< SnvMixOneData > data).counts[1]
    
    def __dealloc__(self):
        free(self._cpt_array)
    
    cdef double * get_resp(self):
        cdef int g
        cdef double * resp
    
        resp = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            resp[g] = self._cpt_array[g]
        
        log_space_normalise_row(resp, NUM_GENOTYPES)
        
        for g in range(NUM_GENOTYPES):
            resp[g] = exp(resp[g])
        
        return resp
    
    cdef double * get_expected_counts_a(self):
        cdef int g
        cdef double * resp, * a
        
        a = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        resp = self.get_resp()
        
        for g in range(NUM_GENOTYPES):
            a[g] = self._a * resp[g]
        
        free(resp)
        
        return a
    
    cdef double * get_expected_counts_b(self):
        cdef int g
        cdef double * resp, * b
        
        b = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        resp = self.get_resp()
        
        for g in range(NUM_GENOTYPES):
            b[g] = self._b * resp[g]
        
        free(resp)
        
        return b
    
    cdef double marginalise(self):
        cdef int g
        cdef double marginal
        
        marginal = 0
        
        for g in range(NUM_GENOTYPES):
            marginal += exp(self._cpt_array[g])
        
        return marginal
    
    cdef _init_cpt_array(self, SnvMixOneData data, SnvMixParameters params):
        cdef int a, b, g
        cdef double mu, log_pi
        
        self._cpt_array = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            a = data.counts[0]
            b = data.counts[1]
            mu = params.mu[g]
            log_pi = log(params.pi[g])
            
            self._cpt_array[g] = log_pi + self._binomial_log_likelihood(a, b, mu)
        
    cdef double _binomial_log_likelihood(self, int a, int b, double mu):
        return a * log(mu) + b * log(1 - mu)
        
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoCpt(SnvMixCpt):
    def __init__(self, SnvMixTwoData data, SnvMixParameters params):
        self._depth = data.depth
        self._make_cpt_array()
        self._init_cpt_array(data, params)
    
    def __dealloc__(self):
        self._free_cpt_array()
    
    cdef double * get_resp(self):
        cdef int g
        cdef double * resp
        
        resp = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):     
            resp[g] = self._resp[g]
        
        return resp
    
    cdef double * get_expected_counts_a(self):
        return self._get_expected_counts(1)
    
    cdef double * get_expected_counts_b(self):
        return self._get_expected_counts(0)
    
    cdef double marginalise(self):
        return self._marginal

    cdef double * _get_expected_counts(self, int a):
        '''
        Get the expected number of counts for match a=1, mismatch a=0
        '''
        cdef int d, g, z
        cdef double read_prob
        cdef double * counts
        
        counts = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            counts[g] = 0
        
        # Marginalise each read over z and d.
        for g in range(NUM_GENOTYPES):
            for d in range(self._depth):
                for z in range(2):
                    counts[g] += self._cpt_array[g][d][a][z]
                                        
        return counts
    
    cdef _init_cpt_array(self, SnvMixTwoData data, SnvMixParameters params):
        self._fill_cpt_array(data, params)
    
    cdef _fill_cpt_array(self, SnvMixTwoData data, SnvMixParameters params):
        cdef int d, g, a, z
        cdef double r, q, mu, pi, read_likelihood, norm_const
        cdef double ** read_marginals, * class_marginals
        
        read_marginals = < double **> malloc(NUM_GENOTYPES * sizeof(double))

        for g in range(NUM_GENOTYPES):
            mu = params.mu[g]
            pi = params.pi[g]
            
            read_marginals[g] = < double *> malloc(self._depth * sizeof(double))
            
            for d in range(self._depth):   
                q = data.q[d]
                r = data.r[d]
                
                read_marginals[g][d] = 0
                
                for a in range(2):
                    for z in range(2):                        
                        read_marginals[g][d] += self._get_read_complete_likelihood(a, z, q, r, mu)
        
        class_marginals = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            class_marginals[g] = params.pi[g]
            
            for d in range(self._depth):
                class_marginals[g] *= read_marginals[g][d]
            
            norm_const += class_marginals[g]
        
        for g in range(NUM_GENOTYPES):
            self._resp[g] = class_marginals[g] / norm_const
    
        for g in range(NUM_GENOTYPES):          
            mu = params.mu[g]
            pi = params.pi[g]

            for d in range(self._depth):
                q = data.q[d]
                r = data.r[d]

                for a in range(2):
                    for z in range(2):                        
                        read_likelihood = self._get_read_complete_likelihood(a, z, q, r, mu)           
                        read_likelihood = read_likelihood * (class_marginals[g] / read_marginals[g][d])
                        self._cpt_array[g][d][a][z] = read_likelihood / norm_const

        self._marginal = norm_const
                    
        for g in range(NUM_GENOTYPES):
            free(read_marginals[g])
        
        free(read_marginals)
        free(class_marginals)                        

    cdef double _get_read_complete_likelihood(self, int a, int z, double q, double r, double mu):
        if a == 0 and z == 0:
            return 0.5 * (1 - r) * (1 - mu)
        elif a == 0 and z == 1:
            return (1 - q) * r * (1 - mu)
        elif a == 1 and z == 0:
            return 0.5 * (1 - r) * mu
        else:
            return q * r * mu
        
    cdef void _make_cpt_array(self):
        cdef int g, d, a, z
        
        self._cpt_array = < double ****> malloc(NUM_GENOTYPES * sizeof(double *))
        
        for g in range(NUM_GENOTYPES):
            self._cpt_array[g] = < double ***> malloc(self._depth * sizeof(double *))
            
            for d in range(self._depth):
                self._cpt_array[g][d] = < double **> malloc(2 * sizeof(double *))
                
                for a in range(2):
                    self._cpt_array[g][d][a] = < double *> malloc(2 * sizeof(double))
                    
                    for z in range(2):
                        self._cpt_array[g][d][a][z] = 0

    cdef void _free_cpt_array(self):
        cdef int d, g, a
        
        for g in range(NUM_GENOTYPES):
            for d in range(self._depth):
                for a in range(2):
                    free(self._cpt_array[g][d][a])
                free(self._cpt_array[g][d])        
            free(self._cpt_array[g])
        free(self._cpt_array)
