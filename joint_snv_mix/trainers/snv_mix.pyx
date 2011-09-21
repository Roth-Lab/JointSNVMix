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
# Models
#=======================================================================================================================       


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
