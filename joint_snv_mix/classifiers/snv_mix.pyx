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
        
        self._mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        self._mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        
        self._pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        self._pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))

    cdef tuple _get_labels(self):
        cdef int normal_genotype, tumour_genotype, joint_genotype 
        cdef list labels
        cdef JointBinaryCounterRow row
        
        row = self._iter._current_row
        
        normal_probs = self._get_probs(row._normal_counts.A, row._normal_counts.B, self._mu_N, self._pi_N)
        tumour_probs = self._get_probs(row._tumour_counts.A, row._tumour_counts.B, self._mu_T, self._pi_T)        
        
        joint_probs = self._combine_probs(normal_probs, tumour_probs)
        
        return joint_probs
                
    cdef tuple _get_probs(self, int a, int b, tuple mu, tuple pi):
        '''
        Return posterior probabilities under SNVMix1 model.
        '''
        cdef int i
        cdef float x
        cdef float probs[3]        
        
        for i in range(3):
            probs[i] = log(pi[i]) + a * log(mu[i]) + b * log(1 - mu[i])
        
        self._normalise_log_probs(probs)
        
        return tuple([x for x in probs[:3]])
        
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
    
    cdef tuple _combine_probs(self, tuple normal_probs, tuple tumour_probs):
        cdef int i, j, k
        cdef float x, total
        cdef float joint_probs[9]
        cdef list normalise_probs
        
        total = 0
        for i in range(3):
            for j in range(3):
                k = 3 * i + j
                joint_probs[k] = normal_probs[i] * tumour_probs[j]
                total += joint_probs[k]
        
        for i in range(9):
            joint_probs[i] = joint_probs[i] / total
        
        return tuple([x for x in joint_probs[:9]])
            
                
        
