'''
Created on 2011-08-04

@author: Andrew Roth
'''
import ConfigParser

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
        free(self.q)
        free(self.r)

cdef SnvMixTwoData makeSnvMixTwoData(base_map_qualities_struct data_struct):
    cdef read_index, i, l
    cdef double temp_q
    cdef SnvMixTwoData data = SnvMixTwoData.__new__(SnvMixTwoData)
    
    l = data_struct.depth.A + data_struct.depth.B
    
    data.depth = l
    
    data.q = < double *> malloc(l * sizeof(double))
    data.r = < double *> malloc(l * sizeof(double))
    
    i = 0
    
    for read_index in range(data_struct.depth.A):
        data.q[i] = get_phred_qual_to_prob(data_struct.base_quals.A[read_index])
        data.r[i] = get_phred_qual_to_prob(data_struct.map_quals.A[read_index])
        
        i += 1
    
    for read_index in range(data_struct.depth.B):        
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
            self._mu[g][0] = mu[g][0]
            self._mu[g][1] = mu[g][1]
            self._pi[g] = pi[g]

    def __str__(self):
        s = "mu_alpha : "
        s += "\t".join([str(x) for x in self.mu_alpha])
        s += "\n"
        
        s += "mu_beta : "
        s += "\t".join([str(x) for x in self.mu_beta])
        s += "\n"
        
        s += "pi : "
        s += "\t".join([str(x) for x in self.pi])
        s += "\n"
        
        return s
    
    property mu_alpha:
        def __get__(self):
            return (self._mu[0][0], self._mu[1][0], self._mu[2][0])
    
    property mu_beta:
        def __get__(self):
            return (self._mu[0][1], self._mu[1][1], self._mu[2][1])
    
    property pi:
        def __get__(self):
            return (self._pi[0], self._pi[1], self._pi[2])

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

    def read_from_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        joint_genotypes = []
        
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        for g in range(NUM_GENOTYPES):
            self._normal_priors._mu[g][0] = float(config.get('mu_N_alpha', genotypes[g]))
            self._tumour_priors._mu[g][0] = float(config.get('mu_T_alpha', genotypes[g]))
            
            self._normal_priors._mu[g][1] = float(config.get('mu_N_beta', genotypes[g]))
            self._tumour_priors._mu[g][1] = float(config.get('mu_T_beta', genotypes[g]))
        
            self._normal_priors._pi[g] = float(config.get('pi_N', genotypes[g]))
            self._tumour_priors._pi[g] = float(config.get('pi_T', genotypes[g]))


#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class SnvMixParameters(object):
    def __init__(self, **kwargs):        
        self._priors = kwargs.get('priors', SnvMixPriors())
        
        mu = kwargs.get('mu', (0.99, 0.5, 0.01))
        pi = kwargs.get('pi', (0.99, 0.009, 0.001))
        
        for g in range(NUM_GENOTYPES):
            self._mu[g] = mu[g]
            self._pi[g] = pi[g]
        
    def __str__(self):
        s = "mu : "
        s += "\t".join([str(x) for x in self.mu]) 
        s += "\n"
        
        s += "pi : "
        s += "\t".join([str(x) for x in self.pi])
        s += "\n"
        
        return s
    
    property mu:
        def __get__(self):
            return tuple([x for x in self._mu[:NUM_GENOTYPES]])
    
    property pi:
        def __get__(self):
            return tuple([x for x in self._pi[:NUM_GENOTYPES]])

    cdef _normalise_pi(self):
        cdef int g
        cdef double norm_const
        
        norm_const = 0
        
        for g in range(NUM_GENOTYPES):
            norm_const += self._pi[g]
        
        for g in range(NUM_GENOTYPES):
            self._pi[g] = self._pi[g] / norm_const

    cdef update(self, double * n, double * a, double * b):
        self._update_mu(a, b)
        self._update_pi(n)

    cdef _update_mu(self, double * a, double * b):
        cdef int g
        cdef double alpha, beta, denom
    
        for g in range(NUM_GENOTYPES):
            alpha = a[g] + self._priors._mu[g][0] - 1
            beta = b[g] + self._priors._mu[g][1] - 1
            denom = alpha + beta

            self._mu[g] = alpha / denom
            
    cdef _update_pi(self, double * n):
        cdef int g
        cdef double denom
        cdef double pi[NUM_GENOTYPES]
        
        denom = 0
        
        for g in range(NUM_GENOTYPES):
            pi[g] = n[g] + self._priors._pi[g] - 1
            denom += pi[g]
        
        for g in range(NUM_GENOTYPES):
            self._pi[g] = pi[g] / denom
    
    cdef double _get_prior_log_likelihood(self):
        cdef double ll
        cdef double x[2]
        
        ll = 0
        
        for g in range(NUM_GENOTYPES):
            x[0] = self._mu[g]
            x[1] = 1 - self._mu[g]
            
            ll += dirichlet_log_likelihood(x, self._priors._mu[g], NUM_BASES)
        
        ll += dirichlet_log_likelihood(self._pi, self._priors._pi, NUM_GENOTYPES)
        
        return ll
            
#---------------------------------------------------------------------------------------------------------------------- 
cdef class PairedSnvMixParameters(object):
    def __init__(self, **kwargs):
        cdef PairedSnvMixPriors priors
        cdef SnvMixPriors normal_priors, tumour_priors
        
        priors = kwargs.get('priors', PairedSnvMixPriors())
        normal_priors = priors._normal_priors
        tumour_priors = priors._tumour_priors
         
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        self._normal_params = SnvMixParameters(priors=normal_priors, mu=mu_N, pi=pi_N)
        self._tumour_params = SnvMixParameters(priors=tumour_priors, mu=mu_T, pi=pi_T)
    
    property normal:
        def __get__(self):
            return self._normal_params
    
    property tumour:
        def __get__(self):
            return self._tumour_params
    
    property mu_N:
        def __get__(self):
            return self._normal_params.mu

    property mu_T:
        def __get__(self):
            return self._tumour_params.mu
        
    property pi_N:
        def __get__(self):
            return self._normal_params.pi

    property pi_T:
        def __get__(self):
            return self._tumour_params.pi

    def write_to_file(self, file_name):
        genotypes = ['AA', 'AB', 'BB']
        
        config = ConfigParser.SafeConfigParser()
        
        config.add_section('pi_N')
        config.add_section('pi_T')
        config.add_section('mu_N')
        config.add_section('mu_T')
        
        for g_N, mu_N in zip(genotypes, self.mu_N):
            config.set('mu_N', g_N, "{0:.10f}".format(mu_N))
        
        for g_T, mu_T in zip(genotypes, self.mu_T):
            config.set('mu_T', g_T, "{0:.10f}".format(mu_T))
        
        for g_N, pi_N in zip(genotypes, self.pi_N):
            config.set('pi_N', g_N, "{0:.10f}".format(pi_N))
        
        for g_T, pi_T in zip(genotypes, self.pi_T):
            config.set('pi_T', g_T, "{0:.10f}".format(pi_T))
        
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
            self._normal_params._mu[g] = float(config.get('mu_N', genotypes[g]))
            self._tumour_params._mu[g] = float(config.get('mu_T', genotypes[g]))
        
            self._normal_params._pi[g] = float(config.get('pi_N', genotypes[g]))
            self._tumour_params._pi[g] = float(config.get('pi_T', genotypes[g]))
        
        self._normal_params._normalise_pi()
        self._tumour_params._normalise_pi()
        
#=======================================================================================================================
# Models
#=======================================================================================================================       
cdef class SnvMixModel(object):
    def __init__(self, SnvMixParameters params):
        self._params = params
    
    def fit(self, list data, convergence_threshold=1e-6, max_iters=100):
        trainer = SnvMixModelTrainer(self, convergence_threshold, max_iters)
        
        trainer.train(data)
        
    property params:
        def __get__(self):
            return self.params

    cdef double _get_lower_bound(self, list data):
        cdef double lb
        cdef SnvMixData pos_data
        
        lb = 0
        
        for pos_data in data:
            lb += self._get_log_likelihood(pos_data)
        
        lb += self._params._get_prior_log_likelihood()
        
        return lb
        
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data):
        pass
    
    cdef double _get_log_likelihood(self, SnvMixData data):
        cdef SnvMixCpt cpt
        cdef double log_likelihood
        
        cpt = self._get_complete_log_likelihood(data)
        
        log_likelihood = cpt.get_log_sum()
        
        return log_likelihood

cdef class PairedSnvMixModel(object):        
    def fit(self, dict data, convergence_threshold=1e-6, max_iters=100):
        self._normal_model.fit(data['normal'], convergence_threshold, max_iters)
        self._tumour_model.fit(data['tumour'], convergence_threshold, max_iters)
    
    property params:
        def __get__(self):
            return self._params

cdef class PairedSnvMixOneModel(PairedSnvMixModel):
    def __init__(self, PairedSnvMixParameters params):
        self._params = params
        
        self._normal_model = SnvMixOneModel(params._normal_params)
        self._tumour_model = SnvMixOneModel(params._tumour_params)

cdef class PairedSnvMixTwoModel(PairedSnvMixModel):
    def __init__(self, PairedSnvMixParameters params):
        self._params = params
        
        self._normal_model = SnvMixTwoModel(params._normal_params)
        self._tumour_model = SnvMixTwoModel(params._tumour_params)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneModel(SnvMixModel):
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data):
        return SnvMixOneCpt(data, self._params)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoModel(SnvMixModel):
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data):
        return SnvMixTwoCpt(data, self._params)
    
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
            print self._model._params

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
        self._model._params.update(ess._n, ess._a, ess._b)

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
            self._a[g] = 0
            self._b[g] = 0
            self._n[g] = 0
    
    cdef update(self, SnvMixCpt cpt):
        cdef double * a, *b, * resp
        
        resp = cpt.get_resp()
        a = cpt.get_expected_counts_a()
        b = cpt.get_expected_counts_b()
    
        for g in range(NUM_GENOTYPES):
            self._n[g] += resp[g]
            self._a[g] += a[g]
            self._b[g] += b[g]
        
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
    
    cdef double get_log_sum(self):
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
    
    cdef double get_log_sum(self):
        cdef int g
        cdef double log_marginal
        
        log_marginal = log_sum_exp(& self._cpt_array[0], NUM_GENOTYPES)
        
        return log_marginal
    
    cdef _init_cpt_array(self, SnvMixOneData data, SnvMixParameters params):
        cdef int a, b, g
        cdef double mu, log_pi
        
        self._cpt_array = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            a = data.counts[0]
            b = data.counts[1]
            mu = params._mu[g]
            log_pi = log(params._pi[g])
            
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
    
    cdef double get_log_sum(self):
        return log(self._marginal)

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
        
        read_marginals = self._get_read_marginals(data, params)
        
        class_marginals = self._get_class_marginals(read_marginals, params)
        
        norm_const = 0
        
        for g in range(NUM_GENOTYPES):
            norm_const += class_marginals[g]
        
        if norm_const == 0:
            norm_const = EPS
        
        for g in range(NUM_GENOTYPES):
            self._resp[g] = class_marginals[g] / norm_const
    
        for g in range(NUM_GENOTYPES):          
            mu = params._mu[g]
            pi = params._pi[g]

            for d in range(self._depth):
                q = data.q[d]
                r = data.r[d]

                for a in range(2):
                    for z in range(2):                        
                        read_likelihood = self._get_read_complete_likelihood(a, z, q, r, mu)           
                        read_likelihood = read_likelihood * (class_marginals[g] / read_marginals[g][d])
                        self._cpt_array[g][d][a][z] = read_likelihood / norm_const

        self._marginal = norm_const
        
        self._free_read_marginals(read_marginals)            
        
        free(class_marginals)                        
    
    cdef double ** _get_read_marginals(self, SnvMixTwoData data, SnvMixParameters params):
        cdef int d, g, a, z
        cdef double ** read_marginals
        
        read_marginals = < double **> malloc(NUM_GENOTYPES * sizeof(double))

        for g in range(NUM_GENOTYPES):
            mu = params._mu[g]
            pi = params._pi[g]
            
            read_marginals[g] = < double *> malloc(self._depth * sizeof(double))
            
            for d in range(self._depth):   
                q = data.q[d]
                r = data.r[d]
                
                read_marginals[g][d] = 0
                
                for a in range(2):
                    for z in range(2):                        
                        read_marginals[g][d] += self._get_read_complete_likelihood(a, z, q, r, mu)
        
        return read_marginals
    
    cdef void _free_read_marginals(self, double ** read_marginals):
        for g in range(NUM_GENOTYPES):
            free(read_marginals[g])
        
        free(read_marginals)
    
    cdef double * _get_class_marginals(self, double ** read_marginals, SnvMixParameters params):
        cdef int d, g, a, z
        cdef double * class_marginals
        
        class_marginals = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            class_marginals[g] = params._pi[g]
            
            for d in range(self._depth):
                class_marginals[g] *= read_marginals[g][d]
        
        return class_marginals

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
