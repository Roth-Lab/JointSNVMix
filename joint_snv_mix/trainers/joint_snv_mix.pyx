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




cdef class JointSnvMixTwoModel(JointSnvMixModel):
    cdef JointSnvMixCpt _get_complete_log_likelihood(self, JointSnvMixData data):
        return JointSnvMixTwoCpt(data, self._params)
    
#=======================================================================================================================
# Trainer
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



