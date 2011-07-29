#cython: cdivision=True

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMixClassifier(Classifier):
    def iter_ref(self, ref, **kwargs):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return SnvMixClassifierRefIterator(
                                           ref,
                                           self._counter.iter_ref(ref),
                                           **kwargs
                                           )
             
cdef class SnvMixClassifierRefIterator(ClassifierRefIterator):
    def __init__(self, char * ref, JointBinaryBaseCounterIterator iter, **kwargs):
        ClassifierRefIterator.__init__(self, ref, iter, **kwargs)
        
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        for i in range(NUM_GENOTYPES):
            self._log_mu_N[i][0] = log(mu_N[i])
            self._log_mu_N[i][1] = log(1 - mu_N[i])
            
            self._log_mu_T[i][0] = log(mu_T[i])
            self._log_mu_T[i][1] = log(1 - mu_T[i])
            
            self._log_pi_N[i] = log(pi_N[i])            
            self._log_pi_T[i] = log(pi_T[i])

    cdef tuple _get_labels(self):
        cdef double x
        cdef double normal_probs[NUM_GENOTYPES]
        cdef double tumour_probs[NUM_GENOTYPES]
        cdef double joint_probs[NUM_JOINT_GENOTYPES] 
        cdef JointBinaryCounterRow row
        
        row = self._iter._current_row
        
        self._get_probs(
                        normal_probs,
                        row._normal_counts.A, row._normal_counts.B,
                        self._log_mu_N, self._log_pi_N
                        )
        self._get_probs(
                        tumour_probs,
                        row._tumour_counts.A, row._tumour_counts.B,
                        self._log_mu_T, self._log_pi_T
                        )        
        
        self._combine_probs(joint_probs, normal_probs, tumour_probs)
        
        return tuple([x for x in joint_probs[:NUM_JOINT_GENOTYPES]])
                                
    cdef void _get_probs(self,
                         double probs[NUM_GENOTYPES],
                         int a,
                         int b,
                         double log_mu[NUM_GENOTYPES][2],
                         double log_pi[NUM_GENOTYPES]):
        '''
        Return posterior probabilities under SNVMix1 model.
        '''
        cdef int i
        cdef double x
        
        for i in range(NUM_GENOTYPES):
            probs[i] = log_pi[i] + a * log_mu[i][0] + b * log_mu[i][1]
        
        log_space_normalise_row(probs, NUM_GENOTYPES)
        
        for i in range(NUM_GENOTYPES):
            probs[i] = exp(probs[i])
    
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
            
                
        
