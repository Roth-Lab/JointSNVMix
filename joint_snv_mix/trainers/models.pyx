cdef class SnvMixModel(object):
    def __init__(self, SnvMixParameters params):
        self._params = params
    
    def fit(self, list data, convergence_threshold=1e-6, max_iters=100):
        trainer = SnvMixModelTrainer(self, convergence_threshold, max_iters)
        
        trainer.train(data)
        
    property params:
        def __get__(self):
            return self.params

    cdef double _get_lower_bound(self, list data):
        cdef double lb
        cdef SnvMixData pos_data
        
        lb = 0
        
        for pos_data in data:
            lb += self._get_log_likelihood(pos_data)
        
        lb += self._params._get_prior_log_likelihood()
        
        return lb
        
    cdef SnvMixCpt _get_complete_log_likelihood(self, SnvMixData data):
        pass
    
    cdef double _get_log_likelihood(self, SnvMixData data):
        cdef SnvMixCpt cpt
        cdef double log_likelihood
        
        cpt = self._get_complete_log_likelihood(data)
        
        log_likelihood = cpt.get_log_sum()
        
        return log_likelihood

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
    
cdef class ModelTrainer(object):    
    def __init__(self, model, convergence_threshold, max_iters):
        self._model = model
        self._convergence_threshold = convergence_threshold
        self._max_iters = max_iters
        
        self._converged = 0
        self._iters = 0
        self._lower_bounds = [FLOAT_INFN]

    cdef train(self, list data):
        while not self._converged:
            ess = self._do_e_step(data)

            self._do_m_step(ess)
            
            self._check_convergence(data)
            
            print self._iters, self._lower_bounds[-1]
            print self._model._params

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
    
cdef class SnvMixModelTrainer(object):    
    cdef SnvMixEss _do_e_step(self, list data):                              
        cdef SnvMixData pos_data
        cdef SnvMixEss ess
        cdef SnvMixCpt cpt
        
        ess = SnvMixEss()
    
        for pos_data in data:
            cpt = self._model._get_complete_log_likelihood(pos_data)                   
            ess.update(cpt)
        
        return ess
    
    cdef void _do_m_step(self, SnvMixEss ess):       
        self._model._params.update(ess._n, ess._a, ess._b)
    
cdef class JointSnvMixModelTrainer(object):
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
