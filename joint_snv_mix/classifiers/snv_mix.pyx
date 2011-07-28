#cython: cdivision=True

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
        
        for i in range(3):
            self._mu_N[i] = mu_N[i]
            self._mu_T[i] = mu_T[i]
            self._pi_N[i] = pi_N[i]
            self._pi_T[i] = pi_T[i]

    cdef tuple _get_labels(self):
        cdef float x
        cdef float normal_probs[3]
        cdef float tumour_probs[3]
        cdef float joint_probs[9] 
        cdef JointBinaryCounterRow row
        
        row = self._iter._current_row
        
        self._get_probs(normal_probs, row._normal_counts.A, row._normal_counts.B, self._mu_N, self._pi_N)
        self._get_probs(tumour_probs, row._tumour_counts.A, row._tumour_counts.B, self._mu_T, self._pi_T)        
        
        self._combine_probs(joint_probs, normal_probs, tumour_probs)
        
        return tuple([x for x in joint_probs[:9]])
                
    cdef void _get_probs(self, float probs[3], int a, int b, float mu[3], float pi[3]):
        '''
        Return posterior probabilities under SNVMix1 model.
        '''
        cdef int i
        cdef float x
        
        for i in range(3):
            probs[i] = log(pi[i]) + a * log(mu[i]) + b * log(1 - mu[i])
        
        self._normalise_log_probs(probs)
        
    cdef void _normalise_log_probs(self, float probs[3]):
        cdef float nc
        
        nc = self._compute_log_norm_constant(probs)
        
        for i in range(3):
            probs[i] = exp(probs[i] - nc)
    
    cdef float _compute_log_norm_constant(self, float probs[3]):
        cdef float max_exp, total
     
        max_exp = probs[0]
     
        for i in range(3):
            if max_exp < probs[i]:
                max_exp = probs[i]
    
        total = 0
        for i in range(3):
            total += exp(probs[i] - max_exp)
        
        return log(total) + max_exp
    
    cdef void _combine_probs(self, float joint_probs[9], float normal_probs[3], float tumour_probs[3]):
        cdef int i, j, k
        cdef float x, total
        cdef list normalise_probs
        
        total = 0
        for i in range(3):
            for j in range(3):
                k = 3 * i + j
                joint_probs[k] = normal_probs[i] * tumour_probs[j]
                total += joint_probs[k]
        
        for i in range(9):
            joint_probs[i] = joint_probs[i] / total
            
                
        
