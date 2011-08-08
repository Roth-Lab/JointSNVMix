'''
Created on 2011-08-04

@author: Andrew Roth
'''
from libc.stdlib cimport free, malloc
from libc.math cimport abs, exp, log

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow
from joint_snv_mix.counters.joint_binary_quality_counter cimport JointBinaryQualityCounterRow
from joint_snv_mix.counters.shared cimport binary_counts_struct, base_map_qualities_struct

from joint_snv_mix.trainers.trainer cimport Trainer
from joint_snv_mix.utils.log_pdf cimport multinomial_log_likelihood, dirichlet_log_likelihood, mixture_posterior

DEF FLOAT_INFN = float('-inf')

DEF NUM_GENOTYPES = 3
DEF NUM_BASES = 2
DEF EPS = 1e-100

#=======================================================================================================================
# Subsamplers
#=======================================================================================================================
cdef class PairedDataSubSampler(object):
    cdef int _skip_size
    cdef int _min_normal_depth
    cdef int _min_tumour_depth
    
    def __init__(self, double skip_size, int min_normal_depth, int min_tumour_depth):
        self._skip_size = skip_size
        self._min_normal_depth = min_normal_depth
        self._min_tumour_depth = min_tumour_depth

    def subsample(self, Counter counter, refs=None):
        cdef int ref_len
        cdef RefIterator ref_iter       
        cdef PairedSampleBinomialCounterRow row
        cdef list sample
        
        sample = {'normal' : [], 'tumour' : []}
        
        print "Randomly sub-sampling every {0}th position the data set.".format(self._skip_size)
        
        if refs = None:
            refs = counter.refs
        
        for ref in refs:
            ref_iter = counter.iter_ref(ref)
            
            i = 0
            
            try:
                while True:
                    ref_iter.cnext()
                    
                    row = ref_iter._current_row
                    
                    if row._normal_depth < self._min_normal_depth or row._tumour_depth < self._min_tumour_depth:
                        continue
                    
                    if i % self._skip_size == 0:               
                        self._add_row_to_sample(sample, row)
                        
                    i += 1
                        
            except StopIteration:
                pass
            
            print "Sub-sampled {0} positions from ref {1}".format(i, ref)
        
        print "Total sub-sample size is {0}".format(len(sample_data['normal']))
        
        return sample_data
    
    cdef _add_row_to_sample(self, dict sample, PairedSampleBinomialCounterRow row):
        pass
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneSubsampler(PairedDataSubSampler):
    cdef _add_row_to_sample(self, dict sample, PairedSampleBinomialCounterRow row):
        cdef SnvMixTwoTrainingData normal_data, tumour_data
    
        normal_data = makeSnvMixTwoTrainingData((< JointBinaryQualityCounterRow > row)._normal_data)
        tumour_data = makeSnvMixTwoTrainingData((< JointBinaryQualityCounterRow > row)._tumour_data)
        
        sample['normal'].append(normal_data)
        sample['tumour'].append(tumour_data)
        
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneSubsampler(PairedDataSubSampler):
    cdef _add_row_to_sample(self, dict sample, PairedSampleBinomialCounterRow row):
        cdef SnvMixOneData normal_data, tumour_data
    
        normal_data = makeSnvMixOneTrainingData((< JointBinaryCounterRow > row)._normal_counts)
        tumour_data = makeSnvMixOneTrainingData((< JointBinaryCounterRow > row)._tumour_counts)
        
        sample['normal'].append(normal_data)
        sample['tumour'].append(tumour_data)

#=======================================================================================================================
# Models
#=======================================================================================================================
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
        
cdef class SnvMixModel(object):
    cdef SnvMixParams params
    
    def __init__(self, SnvMixParams params):
        self.params = params
    
    def fit(self, list data, convergence_threshold=1e-6, max_iters=100):
        trainer = SnvMixModelTrainer(self, convergence_threshold, max_iters)
        
        trainer.train(data)
    
    def predict(self, SnvMixData data):
        pass

    cdef double _get_lower_bound(self, list data):
        cdef double lb
        cdef SnvMixData data
        
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

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneModel(SnvMixModel):
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data)
        return SnvMixOneCpt(data, self.params)
    
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoModel(SnvMixModel):
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data)
        return SnvMixTwoCpt(data, self.params)

#=======================================================================================================================
# CPT
#=======================================================================================================================
cdef class SnvMixCpt(object):
    cdef double * get_resp():
        pass
    
    cdef double * get_expected_counts_a():
        pass
    
    cdef double * get_expected_counts_b():
        pass
    
    cdef double marginalise():
        pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneCpt(SnvMixCpt):
    cdef int _a
    cdef int _b
    cdef double * _cpt_array
    
    def __init__(self, SnvMixData data, SnvMixParams params):
        self._init_cpt_array(data, params)
        
        self._a = (< SnvMixOneData > data).counts[0]
        self._b = (< SnvMixOneData > data).counts[1]
    
    def __dealloc__(self):
        free(self._cpt_array)
    
    cdef double * get_resp():
        cdef int g
        cdef double * resp
    
        resp = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            resp[g] = self._cpt_array[g]
        
        log_space_normalise_row(resp, NUM_GENOTYPES)
        
        return resp
    
    cdef double * get_expected_counts_a()
        cdef int g
        cdef double * resp, * a
        
        a = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        resp = self.get_resp()
        
        for g in range(NUM_GENOTYPES):
            a[g] = self._a * resp[g]
        
        free(resp)
        
        return a
    
    cdef double * get_expected_counts_b()
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
            margianl += exp(self._cpt_array[g])
        
        return marginal
    
    cdef _init_cpt_array(self, SnvMixOneData data, SnvMixParams params):
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
class SnvMixTwoCpt(SnvMixCpt):
    cdef int depth
    cdef double ** _cpt_array
    
    def __init__(self, SnvMixTwoTrainingData data, SnvMixParams params):
        self.depth = data.depth
        
        self._make_cpt_array()
        self._init_cpt_array(SnvMixTwoTrainingData data, SnvMixParams params)
    
    def __dealloc__(self):
        self._destroy_cpt_array()
    
    cdef double * get_resp(self):
        cdef int d, g, a, z
        cdef double read_prob
        cdef double * resp
        cdef double norm_const
        
        resp = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        norm_const = 0
        for g in range(NUM_GENOTYPES):
            resp[g] = 1
        
        # Marginalise each read over a and z.
        for g in range(NUM_GENOTYPES):
            for d in range(self.depth):            
                read_prob = 0       
                for a in range(2):
                    for z in range(2):
                        read_prob += self._cpt_array[d][g][a][z]
                        
                resp[g] = read_prob * resp[g]
            
            norm_const += resp[g]
                
        for g in range(NUM_GENOTYPES):
            resp[g] = resp[g] / norm_const
        
        return resp
    
    cdef double * get_expected_counts_a(self):
        return self._get_expected_counts(1)
    
    cdef double * get_expected_counts_b(self):
        return self._get_expected_counts(0)
    
    cdef double marginalise(self):
        cdef int d, g, a, z
        cdef double marginal, class_marginal
        
        marginal = 0
        
        for g in range(NUM_GENOTYPES):
            class_marginal = 1
            for d in range(self.depth):
                for a in range(2):
                    for z in range(2):
                        class_marginal = class_marginal * self._cpt_array[d][g][a][z]
            
            marginal += class_marginal
        
        return marginal 

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
            for d in range(self.depth):                   
                    for z in range(2):
                        counts[g] += self._cpt_array[d][g][a][z]
        
        return counts
    
    cdef _init_cpt_array(self):
        self._fill_cpt_array()        
        self._normalise_cpt_array()
    
    cdef _fill_cpt_array(self):
        cdef int d, g, a, z
        cdef double r, q, mu, pi, read_likelihood
    
        for d in range(self.depth):
            q = data.q[d]
            r = data.r[d]

            for g in range(NUM_GENOTYPES):            
                mu = params.mu[g]
                pi = params.pi[g]
                
                for a in range(2):
                    for z in range(2):
                        read_likelihood = self._get_read_complete_likelihood(a, z, q, r, mu)                  
                        self._cpt_array[d][g][a][z] = pi * read_likelihood
    
    cdef _normalise_cpt_array(self):
        '''
        Normalise entry in cpt for each read so the sum over a and z adds to 1.
        '''
        cdef int d, g, a, z
        cdef double norm_const
        
        for d in range(self.depth):
            norm_const = 0
            for g in range(NUM_GENOTYPES):
                for a in range(2):
                    for z in range(2):
                        norm_const += self._cpt_array[d][g][a][z]
            
            for g in range(NUM_GENOTYPES):
                for a in range(2):
                    for z in range(2):
                        self._cpt_array[i][g][a][z] = self._cpt_array[d][g][a][z] / norm_const

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
        cdef d, g, a, z
        
        self._cpt_array = < double **> malloc(depth * sizeof(double *))
        
        for d in range(self.depth):
            self._cpt_array[d] = < double **> malloc(NUM_GENOTYPES * sizeof(double *))
            
            for g in range(NUM_GENOTYPES):
                self._cpt_array[d][g] = < double **> malloc(2 * sizeof(double *))
                
                for a in range(2):
                    self._cpt_array[d][g][a] = < double *> malloc(2 * sizeof(double))
                    
                    for z in range(2):
                        self._cpt_array[d][g][a][z] = 0

    cdef void _free_cpt_array(self):
        cdef int d, g, a
        
        for d in range(self.depth):
            for g in range(NUM_GENOTYPES):            
                for a in range(2):
                    free(self._cpt_array[d][g][a])
                free(self._cpt_array[d][g])        
            free(self._cpt_array[d])
        free(self._cpt_array)

#=======================================================================================================================
# Trainer
#=======================================================================================================================
cdef class SnvMixModelTrainer(object):
    cdef bint _converged
    cdef int _iter
    cdef double _convergence_threshold
    cdef int _max_iters

    cdef list _lower_bounds
    cdef SnvMixModel _model
    
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
            params = self._do_m_step(ess)
            
            self._check_convergence(data)
            
            print self._iter, self.lower_bounds[-1]
            print self.model.params

    cdef SnvMixEss _do_e_step(self, list data):                              
        cdef SnvMixTrainingData pos_data
        cdef SnvMixEss ess
        
        ess = SnvMixEss()
    
        for pos_data in data:
            ess.update(pos_data, self._model)
        
        return ess
    
    cdef SnvMixModel _do_m_step(self, SnvMixEss ess):       
        self._model._params._update(ess)

    cdef _check_convergence(self, list data):
        cdef double rel_change, lb
        
        lb = self._model._get_lower_bound(data)        
        self._lower_bounds.append(lb)
        
        rel_change = (lower_bounds[-1] - lower_bounds[-2]) / abs(< double > lower_bounds[-2])
    
        if rel_change < 0:
            print "Lower bound decreased exiting."
            self._converged = 1
        elif rel_change < self._convergence_threshold:
            print "Converged"
            self.converged = 1
        elif self._iter >= self._max_iters:
            print "Maximum number of iters exceeded exiting."
            self._converged = 1
        else:
            self._converged = 0
        
        self._iter += 1

#=======================================================================================================================
# ESS
#=======================================================================================================================
cdef class SnvMixEss(object):
    cdef double a[NUM_GENOTYPES]
    cdef double b[NUM_GENOTYPES]
    cdef double n[NUM_GENOTYPES]
    
    def __init__(self):
        self.reset()
    
    cdef void reset(self):
        cdef int g
        
        for g in range(NUM_GENOTYPES):
            self.a[g] = 0
            self.b[g] = 0
            self.n[g] = 0
    
    cdef update(self, SnvMixData data, SnvMixModel model):
        cdef double * a, *b, * resp
        cdef SnvMixCpt cpt
        
        cpt = model._get_complete_log_likelihood(data)
        
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
# Data
#=======================================================================================================================
cdef class SnvMixData(object):
    pass

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneData(SnvMixData):
    cdef int counts[2]
    
    cdef double component_log_likelihood(self, double mu[2]):
        return multinomial_log_likelihood(self.counts, mu, NUM_BASES)
    
cdef SnvMixOneTrainingData makeSnvMixOneTrainingData(binary_counts_struct counts):
    cdef SnvMixOneTrainingData data = SnvMixOneTrainingData.__new__(SnvMixOneTrainingData)
    
    data.counts[0] = counts.A
    data.counts[1] = counts.B
    
    return data

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoData(SnvMixData):
    cdef int depth
    cdef int * labels
    cdef double * q
    cdef double * r    
    
    def __dealloc__(self):
        free(self.labels)
        free(self.q)
        free(self.r)

cdef SnvMixTwoTrainingData makeSnvMixTwoTrainingData(base_map_qualities_struct data_struct):
    cdef read_index, i, l
    cdef double temp_q
    cdef SnvMixTwoTrainingData data = SnvMixTwoTrainingData.__new__(SnvMixTwoTrainingData)
    
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
# Parameters
#=======================================================================================================================
cdef class PairedSnvMixParameters(object):
    cdef SnvMixParameters _normal_params
    cdef SnvMixParameters _tumour_params

    def __init__(self, **kwargs):
        normal_priors = kwargs.get('priors', PairedSnvMixPriors())
         
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        self._normal_params = SnvMixParameters(priors._normal_priors, mu_N, pi_N)
        self._tumour_params = SnvMixParameters(priors._tumour_priors, mu_T, pi_T)
    
    property normal:
        def __get__(self):
            return self._normal_params
    
    property tumour:
        def __get__(self):
            return self._tumour_params

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixParameters(object):
    cdef SnvMixPriors priors
    
    cdef double mu[NUM_GENOTYPES]
    cdef double pi[NUM_GENOTYPES]

    def __init__(self, SnvMixPriors priors, tuple mu, tuple pi):
        self.priors = priors
        
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

    cdef update(self, SnvMixEss ess):
        self._update_mu(ess)
        self._update_pi(ess)

    cdef _update_mu(self, SnvMixEss ess):
        cdef int g
        cdef double alpha, beta, denom
    
        for g in range(NUM_GENOTYPES)
            alpha = ess.a[g] + self.priors.mu[g][0] - 1
            beta = ess.b[g] + self.priors.mu[g][1] - 1
            denom = alpha + beta
            
            self.mu[g] = alpha / denom
            
    cdef _update_pi(self, SnvMixEss ess):
        cdef int g
        cdef double denom
        cdef double pi[NUM_GENOTYPES]
        
        denom = 0
        
        for g in range(NUM_GENOTYPES):
            pi[g] = ess.n[g] + self.priors.pi[g] - 1
            denom += pi[g]
        
        for g in range(NUM_GENOTYPES):
            self.pi[g] = pi[g] / denom

#=======================================================================================================================
# Priors
#=======================================================================================================================
cdef class PairedSnvMixPriors(object):
    cdef SnvMixPriors _normal_priors
    cdef SnvMixPriors _tumour_priors

    def __init__(self, **kwargs):
        default_mu = (
                      (100, 2),
                      (50, 50),
                      (2, 100)
                      )
        
        default_pi = (1000, 100, 100)
        
        mu_N = kwargs.get('mu_N', default_mu)
        mu_T = kwargs.get('mu_T', default_mu)
                          )
        pi_N = kwargs.get('pi_N', default_pi)
        pi_T = kwargs.get('pi_T', default_pi)
        
        self._normal_priors = SnvMixPriors(mu_N, pi_N)
        self._tumour_priors = SnvMixPriors(mu_T, pi_T)
    
    property normal:
        def __get__(self):
            return self._normal_priors
    
    property tumour:
        def __get__(self):
            return self._tumour_priors

cdef class SnvMixPriors(object):
    cdef double mu[NUM_GENOTYPES][2]
    cdef double pi[NUM_GENOTYPES]

    def __init__(self, tuple mu, tuple pi):
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
