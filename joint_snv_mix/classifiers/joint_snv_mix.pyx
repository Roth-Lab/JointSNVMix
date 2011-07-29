DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class JointSnvMixClassifier(Classifier):
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return JointSnvMixClassifierRefIterator(
                                               ref,
                                               self._counter.iter_ref(ref),
                                               **kwargs
                                               )
             
cdef class JointSnvMixClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryBaseCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        
        pi = kwargs.get('pi', (1e6, 1e2, 1e2, 1e2, 1e4, 1e2, 1e1, 1e1, 1e4))
        
        for i in range(NUM_GENOTYPES):
            self._log_mu_N[i][0] = log(< double > mu_N[i])
            self._log_mu_N[i][1] = log(< double > (1 - mu_N[i]))
            
            self._log_mu_T[i][0] = log(< double > mu_T[i])
            self._log_mu_T[i][1] = log(< double > (1 - mu_T[i]))
        
        nc = sum(pi)
        for i in range(NUM_JOINT_GENOTYPES):            
            self._log_pi[i] = log(< double > pi[i] / nc)

    cdef tuple _get_labels(self):
        cdef double x
        cdef double normal_likelihood[NUM_GENOTYPES]
        cdef double tumour_likelihood[NUM_GENOTYPES]
        cdef double joint_probs[NUM_JOINT_GENOTYPES] 
        cdef JointBinaryCounterRow row
        
        row = self._iter._current_row
        
        self._compute_likelihood(normal_likelihood, row._normal_counts.A, row._normal_counts.B, self._log_mu_N)
        self._compute_likelihood(tumour_likelihood, row._tumour_counts.A, row._tumour_counts.B, self._log_mu_T)   
        
        self._compute_joint_probs(joint_probs, normal_likelihood, tumour_likelihood, self._log_pi)
        
        return tuple([x for x in joint_probs[:NUM_JOINT_GENOTYPES]])
                
    cdef void _compute_likelihood(self, double likelihood[NUM_GENOTYPES], int a, int b, double log_mu[NUM_GENOTYPES][2]):
        '''
        Return posterior probabilities under SNVMix1 model.
        '''
        cdef int i
        
        for i in range(NUM_GENOTYPES):
            likelihood[i] = a * log_mu[i][0] + b * log_mu[i][1]
    
    cdef void _compute_joint_probs(self,
                                   double joint_probs[NUM_JOINT_GENOTYPES],
                                   double normal_likelihood[NUM_GENOTYPES],
                                   double tumour_likelihood[NUM_GENOTYPES],
                                   double log_pi[NUM_JOINT_GENOTYPES]):
        cdef int i, j, k
        cdef double x, total
        
        total = 0
        for i in range(NUM_GENOTYPES):
            for j in range(NUM_GENOTYPES):
                k = NUM_GENOTYPES * i + j
                joint_probs[k] = log_pi[i] + normal_likelihood[i] + tumour_likelihood[j]

        log_space_normalise_row(joint_probs, NUM_JOINT_GENOTYPES)
        
        for i in range(NUM_JOINT_GENOTYPES):
            joint_probs[i] = exp(joint_probs[i])
