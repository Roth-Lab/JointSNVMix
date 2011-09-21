'''
Created on 2011-08-04

@author: Andrew Roth
'''
import ConfigParser

DEF FLOAT_INFN = float('-inf')

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9
DEF NUM_BASES = 2
DEF EPS = 1e-100

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
        cdef double log_likelihood
        
        cpt = self._get_complete_log_likelihood(data)
        
        log_likelihood = cpt.get_log_sum()
        
        return log_likelihood

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
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._n[g] += resp[g]
        
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
    
    cdef double get_log_sum(self):
        pass

#---------------------------------------------------------------------------------------------------------------------- 
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
        
#---------------------------------------------------------------------------------------------------------------------- 
cdef class JointSnvMixTwoCpt(JointSnvMixCpt):
    def __init__(self, SnvMixTwoData normal_data, SnvMixTwoData tumour_data, JointSnvMixParameters params):
        cdef SampleCpt normal_cpt, tumour_cpt
        cdef double * joint_class_marginals
        
        normal_cpt = makeSampleCpt(normal_data, params._mu_N)
        tumour_cpt = makeSampleCpt(tumour_data, params._mu_T)
        
        joint_class_marginals = self._get_joint_class_marginals(normal_cpt._class_marginals,
                                                                tumour_cpt._class_marginals,
                                                                params._pi)
                
        self._init_marginal(joint_class_marginals)
        self._init_resp(joint_class_marginals)
        
        self._init_normal_expected_counts(normal_cpt, joint_class_marginals)
        self._init_tumour_expected_counts(tumour_cpt, joint_class_marginals)
        
        free(joint_class_marginals)
    
    def __dealloc__(self):
        free(self._normal_counts_a)
        free(self._normal_counts_b)
        free(self._tumour_counts_a)
        free(self._tumour_counts_b)

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
    
    cdef double get_log_sum(self):
        return log(self._marginal)
    
    cdef double * _get_joint_class_marginals(self, double * normal_marginals, double * tumour_marginals, double * pi):
        cdef int g_N, g_T, g_J
        cdef double * joint_marginals
        
        joint_marginals = < double *> malloc(NUM_JOINT_GENOTYPES * sizeof(double))
        
        for g_N in range(NUM_GENOTYPES):
            for g_T in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                
                joint_marginals[g_J] = pi[g_J] * normal_marginals[g_N] * tumour_marginals[g_T]
    
        return joint_marginals
    
    cdef void _init_marginal(self, double * joint_marginals):
        cdef int g_J
        
        self._marginal = 0
        
        for g_J in range(NUM_JOINT_GENOTYPES):
            self._marginal += joint_marginals[g_J]
        
        if self._marginal == 0:
            self._marginal = EPS
    
    cdef void _init_resp(self, double * joint_marginals):
        cdef int g_J
        
        for g_J in range(NUM_JOINT_GENOTYPES):                
                self._resp[g_J] = joint_marginals[g_J] / self._marginal
    
    cdef void _init_normal_expected_counts(self, SampleCpt cpt, double * joint_marginals):
        cdef int g_N, g_T, d
        cdef double norm_const[NUM_GENOTYPES]
        
        for g_N in range(NUM_GENOTYPES):
            norm_const[g_N] = 0
            for g_T in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                        
                norm_const[g_N] += joint_marginals[g_J]
            
            norm_const[g_N] = norm_const[g_N] / self._marginal 
        
        self._normal_counts_a = cpt.get_expected_counts_a(norm_const)
        self._normal_counts_b = cpt.get_expected_counts_b(norm_const)

    cdef void _init_tumour_expected_counts(self, SampleCpt cpt, double * joint_marginals):
        cdef int g_N, g_T, d
        cdef double norm_const[NUM_GENOTYPES]
        
        for g_T in range(NUM_GENOTYPES):
            norm_const[g_T] = 0
            for g_N in range(NUM_GENOTYPES):
                g_J = NUM_GENOTYPES * g_N + g_T
                        
                norm_const[g_T] += joint_marginals[g_J]
            
            norm_const[g_T] = norm_const[g_T] / self._marginal 

        self._tumour_counts_a = cpt.get_expected_counts_a(norm_const)
        self._tumour_counts_b = cpt.get_expected_counts_b(norm_const)

cdef SampleCpt makeSampleCpt(SnvMixTwoData data, double * mu):
    cdef SampleCpt cpt = SampleCpt.__new__(SampleCpt)
    
    cpt._depth = data.depth
    
    cpt._init_cpt_array(data, mu)
    cpt._init_read_marginals()
    cpt._init_class_marginals()
    
    return cpt
    

cdef class SampleCpt(object):
    def __dealloc__(self):
        self._free_cpt_array()
        self._free_read_marginals()
        free(self._class_marginals)

    cdef double * get_expected_counts_a(self, double * norm_const):
        return self._get_expected_counts(1, norm_const)

    cdef double * get_expected_counts_b(self, double * norm_const):
        return self._get_expected_counts(0, norm_const)
    
    cdef double * _get_expected_counts(self, int a, double * norm_const):
        cdef int g, d
        cdef double read_prob
        cdef double * counts
        
        counts = < double *> malloc(NUM_GENOTYPES * sizeof(double))
    
        for g in range(NUM_GENOTYPES):
            counts[g] = 0
            
            for d in range(self._depth):   
                read_prob = norm_const[g] / self._read_marginals[g][d]
                counts[g] += read_prob * (self._cpt_array[g][d][a][0] + self._cpt_array[g][d][a][1])
        
        return counts

    cdef void _init_cpt_array(self, SnvMixTwoData data, double * mu):
        cdef int g, d, a, z
        cdef double m, q, r 

        self._allocate_cpt_array()

        for g in range(NUM_GENOTYPES):
            m = mu[g]
            
            for d in range(self._depth):   
                q = data.q[d]
                r = data.r[d]
                
                for a in range(2):
                    for z in range(2):                        
                        self._cpt_array[g][d][a][z] = self._get_read_complete_likelihood(a, z, q, r, m)

    cdef double _get_read_complete_likelihood(self, int a, int z, double q, double r, double mu):
        if a == 0 and z == 0:
            return 0.5 * (1 - r) * (1 - mu)
        elif a == 0 and z == 1:
            return (1 - q) * r * (1 - mu)
        elif a == 1 and z == 0:
            return 0.5 * (1 - r) * mu
        else:
            return q * r * mu                          
    
    cdef void _init_read_marginals(self):
        cdef int  g, d, a, z
        
        self._allocate_read_marginals()

        for g in range(NUM_GENOTYPES):
            for d in range(self._depth):
                self._read_marginals[g][d] = 0
                
                for a in range(2):
                    for z in range(2):                        
                        self._read_marginals[g][d] += self._cpt_array[g][d][a][z]
    
    cdef void _init_class_marginals(self):
        cdef int  g, d, a, z
        
        self._class_marginals = < double *> malloc(NUM_GENOTYPES * sizeof(double))
        
        for g in range(NUM_GENOTYPES):
            self._class_marginals[g] = 1
            
            for d in range(self._depth):
                self._class_marginals[g] * = self._read_marginals[g][d]
    
    cdef void _allocate_cpt_array(self):
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
    
    cdef void _allocate_read_marginals(self):
        cdef int g
         
        self._read_marginals = < double **> malloc(NUM_GENOTYPES * sizeof(double))

        for g in range(NUM_GENOTYPES):
            self._read_marginals[g] = < double *> malloc(self._depth * sizeof(double))
        

    cdef void _free_read_marginals(self):
        for g in range(NUM_GENOTYPES):
            free(self._read_marginals[g])
        
        free(self._read_marginals)

