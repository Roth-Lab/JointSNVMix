cdef class JointSnvMixOneModel(JointSnvMixModel):
    cdef JointSnvMixCpt _get_complete_log_likelihood(self, JointSnvMixData data):
        return JointSnvMixOneCpt(data, self._params)

cdef class JointSnvMixOneCpt(JointSnvMixCpt):
    def __init__(self, SnvMixOneData normal_data, SnvMixOneData tumour_data, JointSnvMixParameters params):
        self._init_cpt_array(data, params)
        
        self._a_N = normal_data.counts[0]
        self._b_N = normal_data.counts[1]
        
        self._a_T = tumour_data.counts[0]
        self._b_T = tumour_data.counts[1]

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
        
        marginal = 0
        
        log_marginal = log_sum_exp(& self._cpt_array[0], NUM_JOINT_GENOTYPES)        
        
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