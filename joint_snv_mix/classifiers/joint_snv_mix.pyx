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
        
        for i in range(3):
            self._log_mu_N[i][0] = log(< float > mu_N[i])
            self._log_mu_N[i][1] = log(< float > (1 - mu_N[i]))
            
            self._log_mu_T[i][0] = log(< float > mu_T[i])
            self._log_mu_T[i][1] = log(< float > (1 - mu_T[i]))
        
        nc = sum(pi)
        for i in range(9):            
            self._log_pi[i] = log(< float > pi[i] / nc)

    cdef tuple _get_labels(self):
        cdef float x
        cdef float normal_likelihood[3]
        cdef float tumour_likelihood[3]
        cdef float joint_probs[9] 
        cdef JointBinaryCounterRow row
        
        row = self._iter._current_row
        
        self._compute_likelihood(normal_likelihood, row._normal_counts.A, row._normal_counts.B, self._log_mu_N)
        self._compute_likelihood(tumour_likelihood, row._tumour_counts.A, row._tumour_counts.B, self._log_mu_T)   
        
        self._compute_joint_probs(joint_probs, normal_likelihood, tumour_likelihood, self._log_pi)
        
        return tuple([x for x in joint_probs[:9]])
                
    cdef void _compute_likelihood(self, float likelihood[3], int a, int b, float log_mu[3][2]):
        '''
        Return posterior probabilities under SNVMix1 model.
        '''
        cdef int i
        
        for i in range(3):
            likelihood[i] = a * log_mu[i][0] + b * log_mu[i][1]
    
    cdef void _compute_joint_probs(self,
                                   float joint_probs[9],
                                   float normal_likelihood[3],
                                   float tumour_likelihood[3],
                                   float log_pi[9]):
        cdef int i, j, k
        cdef float x, total
        
        total = 0
        for i in range(3):
            for j in range(3):
                k = 3 * i + j
                joint_probs[k] = log_pi[i] + normal_likelihood[i] + tumour_likelihood[j]

        self._normalise_log_probs(joint_probs)

    cdef void _normalise_log_probs(self, float probs[9]):
        cdef float nc
        
        nc = self._compute_log_norm_constant(probs)
        
        for i in range(9):
            probs[i] = exp(probs[i] - nc)
    
    cdef float _compute_log_norm_constant(self, float probs[9]):
        cdef float max_exp, total
     
        max_exp = probs[0]
     
        for i in range(9):
            if max_exp < probs[i]:
                max_exp = probs[i]
    
        total = 0
        for i in range(9):
            total += exp(probs[i] - max_exp)
        
        return log(total) + max_exp
