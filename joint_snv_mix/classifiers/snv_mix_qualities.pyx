#cython: cdivision=True

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMix2Classifier(Classifier):
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return SnvMix2ClassifierRefIterator(
                                           ref,
                                           self._counter.iter_ref(ref),
                                           **kwargs
                                           )
             
cdef class SnvMix2ClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryQualityCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        self._init_qual_map()
        
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        for i in range(NUM_GENOTYPES):
            self._mu_N[i] = mu_N[i]            
            self._mu_T[i] = mu_T[i]
            
            self._log_pi_N[i] = log(pi_N[i])            
            self._log_pi_T[i] = log(pi_T[i])
    
    cdef _init_qual_map(self):
        cdef int qual
        cdef double base, exp
        
        for qual in range(256):
            exp = -1 * (< double > qual) / 10
            base = 10
            
            self._qual_map[qual] = 1 - pow(base, exp)

    cdef tuple _get_labels(self):
        cdef double x
        cdef double normal_probs[NUM_GENOTYPES]
        cdef double tumour_probs[NUM_GENOTYPES]
        cdef double joint_probs[NUM_JOINT_GENOTYPES] 
        cdef JointBinaryQualityCounterRow row
        
        row = self._iter._current_row
        
        self._get_probs(
                        normal_probs,
                        row._normal_data,
                        self._mu_N, self._log_pi_N
                        )
        self._get_probs(
                        tumour_probs,
                        row._tumour_data,
                        self._mu_T, self._log_pi_T
                        )        
        
        self._combine_probs(joint_probs, normal_probs, tumour_probs)
        
        return tuple([x for x in joint_probs[:NUM_JOINT_GENOTYPES]])
                                
    cdef void _get_probs(self,
                         double probs[NUM_GENOTYPES],
                         base_map_qualities_struct data,
                         double mu[NUM_GENOTYPES],
                         double log_pi[NUM_GENOTYPES]):
        '''
        Return posterior probabilities under SNVMix2 model.
        '''
        cdef int i, genotype
        cdef double p, q, r, m, aligned_wrong, aligned_right
        
        for genotype in range(NUM_GENOTYPES):
            p = 0            
            m = mu[genotype]
            
            for read_index in range(data.depth.A):
                q = self._qual_map[data.base_quals.A[read_index]]
                r = self._qual_map[data.map_quals.A[read_index]]
                
                p += self._compute_single_base_log_prob(q, r, m)
            
            for read_index in range(data.depth.B):
                q = (1 - self._qual_map[data.base_quals.B[read_index]]) / 3
                r = self._qual_map[data.map_quals.B[read_index]]
                
                p += self._compute_single_base_log_prob(q, r, m)
            
            probs[genotype] = log_pi[genotype] + p 
        
        log_space_normalise_row(probs, NUM_GENOTYPES)
        
        for i in range(NUM_GENOTYPES):
            probs[i] = exp(probs[i]) 
            
    cdef double _compute_single_base_log_prob(self, double q, double r, double m):
        aligned_wrong = 0.5 * (1 - r)
        aligned_right = r * ((1 - q) * (1 - m) + q * m)
        
        return log(aligned_wrong + aligned_right)
    
    cdef void _combine_probs(self,
                             double joint_probs[NUM_JOINT_GENOTYPES],
                             double normal_probs[NUM_GENOTYPES],
                             double tumour_probs[NUM_GENOTYPES]):
        cdef int i, j, k
        cdef double total
        
        total = 0
        for i in range(NUM_GENOTYPES):
            for j in range(NUM_GENOTYPES):
                k = NUM_GENOTYPES * i + j
                
                joint_probs[k] = normal_probs[i] * tumour_probs[j]

                total += joint_probs[k]
                
        for i in range(NUM_JOINT_GENOTYPES):
            joint_probs[i] = joint_probs[i] / total         
        
