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