'''
Created on 2011-08-04

@author: Andrew Roth
'''
from libc.stdlib cimport free, malloc
from libc.math cimport exp, log

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.joint_binary_counter cimport JointBinaryCounterRow

from joint_snv_mix.trainers.trainer cimport Trainer
from joint_snv_mix.utils.log_pdf cimport multinomial_log_likelihood, dirichlet_log_likelihood, mixture_posterior

DEF FLOAT_INFN = float('-inf')

DEF NUM_GENOTYPES = 3
DEF NUM_BASES = 2

ctypedef struct snv_mix_params_struct:
    double mu[NUM_GENOTYPES][NUM_BASES]
    double pi[NUM_GENOTYPES]

ctypedef struct snv_mix_priors_struct:
    double mu[NUM_GENOTYPES][NUM_BASES]
    double pi[NUM_GENOTYPES]
    
ctypedef struct snv_mix_ess_struct:
    double n[NUM_GENOTYPES]
    double a[NUM_GENOTYPES]
    double b[NUM_GENOTYPES]
    
cdef class SnvMixParameters(object):
    cdef snv_mix_params_struct _normal
    cdef snv_mix_params_struct _tumour

    def __init__(self, **kwargs):
        self._init_params(**kwargs)
    
    def _init_params(self, **kwargs):
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        for i in range(NUM_GENOTYPES):
            self._normal.mu[i][0] = mu_N[i]
            self._normal.mu[i][1] = 1 - mu_N[i]
            
            self._tumour.mu[i][0] = mu_T[i]
            self._tumour.mu[i][1] = 1 - mu_T[i]
            
            self._normal.pi[i] = pi_N[i]            
            self._tumour.pi[i] = pi_T[i]
    
    property normal:
        def __get__(self):
            return self._convert_params_struct_to_dict(self._normal)
    
    property tumour:
        def __get__(self):
            return self._convert_params_struct_to_dict(self._tumour)
    
    cdef dict _convert_params_struct_to_dict(self, snv_mix_params_struct params):
        result = {}
        
        result['pi'] = []
        result['mu'] = []
        
        for i in range(NUM_GENOTYPES):
            result['pi'].append(params.pi[i])
            result['mu'].append(params.mu[i][0])
        
        return result

cdef SnvMixParameters makeSnvMixParameters(snv_mix_params_struct normal_params, snv_mix_params_struct tumour_params):
     cdef SnvMixParameters params = SnvMixParameters.__new__(SnvMixParameters)
    
     params._normal = normal_params
     params._tumour = tumour_params
     
     return params

cdef class SnvMixPriors(object):
    cdef snv_mix_priors_struct _normal
    cdef snv_mix_priors_struct _tumour

    def __init__(self, **kwargs):
        self._init_params(**kwargs)
    
    def _init_params(self, **kwargs):
        mu_N = kwargs.get('mu_N', (
                                   (100, 2),
                                   (50, 50),
                                   (2, 100)
                                   )
                          )
        mu_T = kwargs.get('mu_T', (
                                   (100, 2),
                                   (50, 50),
                                   (2, 100)
                                   )
                          )
        pi_N = kwargs.get('pi_N', (1000, 100, 100))
        pi_T = kwargs.get('pi_T', (1000, 100, 100))
        
        for i in range(NUM_GENOTYPES):
            self._normal.mu[i][0] = mu_N[i][0]
            self._normal.mu[i][1] = mu_N[i][1]
            
            self._tumour.mu[i][0] = mu_T[i][0]
            self._tumour.mu[i][1] = mu_T[i][1]
            
            self._normal.pi[i] = pi_N[i]            
            self._tumour.pi[i] = pi_T[i]
    
    property normal:
        def __get__(self):
            return self._convert_priors_struct_to_dict(self._normal)
    
    property tumour:
        def __get__(self):
            return self._convert_priors_struct_to_dict(self._tumour)
    
    cdef dict _convert_priors_struct_to_dict(self, snv_mix_priors_struct params):
        result = {}
        
        result['pi'] = []
        result['mu'] = []
        
        for i in range(NUM_GENOTYPES):
            result['pi'].append(params.pi[i])
            result['mu'].append(
                                (params.mu[i][0], params.mu[i][1])
                                )
        
        return result
    

cdef class SnvMixTrainer(Trainer):
    def train(self, Counter counter, init_params, priors, sub_sample_fraction=0.01):
        sample_data = self._load_data_set(counter, sub_sample_fraction)
        
        print "training"
        
        parameters = self._train(sample_data, init_params, priors)
        
        return parameters
    
    cdef _train(self, list sample_data, SnvMixParameters init_params, SnvMixPriors priors):
        cdef int i, l
        cdef int ** normal_counts
        cdef int ** tumour_counts
        cdef snv_mix_params_struct normal_params, tumour_params
        cdef JointBinaryCounterRow row
        
        l = len(sample_data)

        normal_counts = < int **> malloc(l * sizeof(int *))
        tumour_counts = < int **> malloc(l * sizeof(int *))
        
        for i, row in enumerate(sample_data):
            normal_counts[i] = < int *> malloc(2 * sizeof(int))
            normal_counts[i][0] = row._normal_counts.A
            normal_counts[i][1] = row._normal_counts.B
            
            tumour_counts[i] = < int *> malloc(2 * sizeof(int))
            tumour_counts[i][0] = row._tumour_counts.A
            tumour_counts[i][1] = row._tumour_counts.B     
    
        normal_params = self._do_em(normal_counts, init_params._normal, priors._normal, l)
        tumour_params = self._do_em(tumour_counts, init_params._tumour, priors._tumour, l)
        
        for i in range(l):
            free(normal_counts[i])
            free(tumour_counts[i])
        
        free(normal_counts)
        free(tumour_counts)
        
        return makeSnvMixParameters(normal_params, tumour_params)
    
    cdef snv_mix_params_struct _do_em(self,
                                      int ** counts,
                                      snv_mix_params_struct init_params,
                                      snv_mix_priors_struct priors,
                                      int sample_size):
        cdef bint converged
        cdef int iter
        cdef double lb
        cdef list lower_bounds
        cdef snv_mix_params_struct params
        cdef snv_mix_ess_struct ess
        
        iter = 0
        converged = 0
        lower_bounds = [FLOAT_INFN]
        
        params = init_params
        
        lb = self._get_lower_bound(counts, params, priors, sample_size)        
        lower_bounds.append(lb)
                
        converged = self._check_convergence(lower_bounds, iter)
        
        while not converged:
            ess = self._do_e_step(counts, params, sample_size)
            params = self._do_m_step(ess, priors)
            
            lb = self._get_lower_bound(counts, params, priors, sample_size)        
            lower_bounds.append(lb)
            
            iter += 1
            
            converged = self._check_convergence(lower_bounds, iter)
            
            print iter, lb
            print "mu : ", params.mu[0][0], params.mu[1][0], params.mu[2][0]
            print "pi : ", params.pi[0], params.pi[1], params.pi[2]
        
        return params
    
    cdef snv_mix_ess_struct _do_e_step(self,
                    (int *) counts[2],
                    snv_mix_params_struct params,
                    int sample_size):
        cdef int i, g
        cdef double ll[NUM_GENOTYPES], reps[NUM_GENOTYPES]
        cdef snv_mix_ess_struct ess
        
        ess = make_snv_mix_ess_struct()
    
        for i in range(sample_size):
            for g in range(NUM_GENOTYPES):
                ll[g] = multinomial_log_likelihood(counts[i], params.mu[g], NUM_BASES)
            
            resp = mixture_posterior(ll, params.pi, NUM_BASES)            
            
            for g in range(NUM_GENOTYPES):
                ess.n[g] += resp[g]
                ess.a[g] += counts[i][0] * resp[g]
                ess.b[g] += counts[i][1] * resp[g]
        
        return ess
    
    cdef snv_mix_params_struct _do_m_step(self,
                                          snv_mix_ess_struct ess,
                                          snv_mix_priors_struct priors):
        cdef int g
        cdef double alpha, beta, pi, denom
        cdef snv_mix_params_struct params
    
        for g in range(NUM_GENOTYPES):
            alpha = ess.a[g] + priors.mu[g][0] - 1
            beta = ess.b[g] + priors.mu[g][1] - 1
            denom = alpha + beta
            
            params.mu[g][0] = alpha / denom
            params.mu[g][1] = 1 - params.mu[g][0]
        
        denom = 0
        for g in range(NUM_GENOTYPES):
            pi = ess.n[g] + priors.pi[g] - 1
            denom += pi
        
        for g in range(NUM_GENOTYPES):
            pi = ess.n[g] + priors.pi[g] - 1
            
            params.pi[g] = pi / denom
        
        return params
    
    cdef double _get_lower_bound(self,
                                 int ** counts,
                                 snv_mix_params_struct params,
                                 snv_mix_priors_struct priors,
                                 int sample_size):
        cdef int i, g
        cdef double ll, pos_likelihod
        
        ll = 0
        
        for i in range(sample_size):
            pos_likelihood = 0
            for g in range(NUM_GENOTYPES):
                pos_likelihood += exp(log(params.pi[g]) + multinomial_log_likelihood(counts[i], params.mu[g], NUM_BASES))
            
            ll += log(pos_likelihood)
        
        for g in range(NUM_GENOTYPES):
            ll += dirichlet_log_likelihood(params.mu[g], priors.mu[g], NUM_BASES)
        
        ll += dirichlet_log_likelihood(params.pi, priors.pi, NUM_GENOTYPES)
        
        return ll
    
    cdef bint _check_convergence(self, lower_bounds, iter):
        if lower_bounds[-1] - lower_bounds[-2] < 0:
            print "Lower bound decreased exiting."
            return 1
        elif lower_bounds[-1] - lower_bounds[-2] < 1e-6:
            print "Converged"
            return 1
        elif iter >= 100:
            print "Maximum number of iters exceeded exiting."
            return 1
        else:
            return 0

#=======================================================================================================================
# Helper functions
#=======================================================================================================================
cdef snv_mix_ess_struct make_snv_mix_ess_struct():
    cdef int g
    cdef snv_mix_ess_struct ess
    
    for g in range(NUM_GENOTYPES):
        ess.n[g] = 0
        ess.a[g] = 0
        ess.b[g] = 0
    
    return ess
