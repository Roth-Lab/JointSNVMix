NUM_GENOTYPES = 3

cdef class Parameters(object):
    def __str__(self):
        for param_name in self.parameters:
            print param_name, " : ",
            print "\t".join(str([x for x in getattr(self, params_name)]))

    def read_from_file(self, file_name):
        parameters = {}
        
        config = ConfigParser.SafeConfigParser()
        config.read(file_name)
        
        for param_name in self.__dict__:
            values = []
                        
            for genotype in enumerate(sections[section]):
                param = float(config.get(sections, genotype))
                
                parameters[section][genotype] = param
        
        self.parameters = parameters 
    
    def write_to_file(self, file_name):
        config = ConfigParser.SafeConfigParser()
        
        for param_name in self.parameters:
            for genotype in self.parameters[param_name]:
                value = self.parameters[param_name][genotype]
                
                config.set(param_name, genotype, "{0:.10f}".format(value))
        
        fh = open(file_name, 'w')
        config.write(fh)
        fh.close()    

#=======================================================================================================================
# Priors
#=======================================================================================================================
cdef class JointSnvMixHyperParameters(Parameters):
    def __init__(self, **kwargs):
        default_mu = (
                      (100, 2),
                      (50, 50),
                      (2, 100)
                      )
        
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        
        self.mu_N = kwargs.get('mu_N', default_mu)
        self.mu_T = kwargs.get('mu_T', default_mu)
        self.pi = kwargs.get('pi', default_pi)

    propert mu_N:
        def __get__(self):
            return tuple([(self._mu_N[x][0], self._mu_N[x][1]) for x in range(NUM_GENOTYPES)]])
        
        def __set__(self, value):
            for g in range(NUM_GENOTYPES):
                self._mu_N[g][0] = value[g][0]
                self._mu_N[g][1] = value[g][1]

    propert mu_T:
        def __get__(self):
            return tuple([(self._mu_T[x][0], self._mu_T[x][1]) for x in range(NUM_GENOTYPES)]])
        
        def __set__(self, value):
            for g in range(NUM_GENOTYPES):
                self._mu_T[g][0] = value[g][0]
                self._mu_T[g][1] = value[g][1]

    propert pi:
        def __get__(self):
            return tuple([x for x in self._pi[:NUM_JOINT_GENOTYPES]])
        
        def __set__(self, value):
            for g in range(NUM_JOINT_GENOTYPES):
                self._pi[g] = value[g]

#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class JointSnvMixParameters(Parameters):
    def __init__(self, **kwargs):        
        self._priors = kwargs.get('priors', JointSnvMixPriors())
        
        default_mu = (0.99, 0.5, 0.01)
        default_pi = (1e6, 1e3, 1e3, 1e3, 1e4, 1e3, 1e1, 1e1, 1e4)
        
        mu_N = kwargs.get('mu_N', default_mu)
        mu_T = kwargs.get('mu_T', default_mu)        
        pi = kwargs.get('pi', default_pi) 
            
    property mu_N:
        def __get__(self):
            return tuple([x for x in self._mu_N[:NUM_GENOTYPES]])
        
        def __set__(self, value):
            for g in range(NUM_GENOTYPES):
                self._mu_N[g][0] = value[g]   
    
    property mu_T:
        def __get__(self):
            return tuple([x for x in self._mu_T[:NUM_GENOTYPES]])
        
        def __set__(self, value):
            for g in range(NUM_GENOTYPES):
                self._mu_T[g][0] = value[g]           
    
    property pi:
        def __get__(self):
            return tuple([x for x in self._pi[:NUM_JOINT_GENOTYPES]])
        
        def __set__(self, value):
            for g in range(NUM_JOINT_GENOTYPES):
                self._pi[g] = value[g]
            
            self._normalise_pi()
    
    cdef _normalise_pi(self):
        cdef int g
        cdef double norm_const
        
        norm_const = 0
        
        for g in range(NUM_JOINT_GENOTYPES):
            norm_const += self._pi[g]
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = self._pi[g] / norm_const

    cdef update(self, double * n, double * a_N, double * a_T, double * b_N, double * b_T):
        self._update_mu(self._mu_N, self._priors._mu_N, a_N, b_N)
        self._update_mu(self._mu_T, self._priors._mu_T, a_T, b_T)
        self._update_pi(n)

    cdef _update_mu(self, double * mu, double mu_prior[NUM_GENOTYPES][2], double * a, double * b):
        cdef int g
        cdef double alpha, beta, denom
    
        for g in range(NUM_GENOTYPES):
            alpha = a[g] + mu_prior[g][0] - 1
            beta = b[g] + mu_prior[g][1] - 1
            denom = alpha + beta

            mu[g] = alpha / denom
            
    cdef _update_pi(self, double * n):
        cdef int g
        cdef double denom
        cdef double pi[NUM_JOINT_GENOTYPES]
        
        denom = 0
        
        for g in range(NUM_JOINT_GENOTYPES):
            pi[g] = n[g] + self._priors._pi[g] - 1
            denom += pi[g]
        
        for g in range(NUM_JOINT_GENOTYPES):
            self._pi[g] = pi[g] / denom
    
    cdef double _get_prior_log_likelihood(self):
        cdef double ll
        cdef double x[2]
        
        ll = 0
        
        for g in range(NUM_GENOTYPES):
            x[0] = self._mu_N[g]
            x[1] = 1 - self._mu_N[g]            
            ll += dirichlet_log_likelihood(x, self._priors._mu_N[g], NUM_BASES)
            
            x[0] = self._mu_T[g]
            x[1] = 1 - self._mu_T[g]
            ll += dirichlet_log_likelihood(x, self._priors._mu_T[g], NUM_BASES)
        
        ll += dirichlet_log_likelihood(self._pi, self._priors._pi, NUM_JOINT_GENOTYPES)
        
        return ll
    
#=======================================================================================================================
# Priors
#=======================================================================================================================
cdef class SnvMixHyperParameters(Parameters):
    def __init__(self, **kwargs):
        default_mu = (
                      (100, 2),
                      (50, 50),
                      (2, 100)
                      )
        
        default_pi = (1000, 100, 100)
        
        mu = kwargs.get('mu', default_mu)
        pi = kwargs.get('pi', default_pi)
        
        for g in range(NUM_GENOTYPES):
            self._mu[g][0] = mu[g][0]
            self._mu[g][1] = mu[g][1]
            self._pi[g] = pi[g]
    
    property mu_alpha:
        def __get__(self):
            return (self._mu[0][0], self._mu[1][0], self._mu[2][0])
    
    property mu_beta:
        def __get__(self):
            return (self._mu[0][1], self._mu[1][1], self._mu[2][1])
    
    property pi:
        def __get__(self):
            return (self._pi[0], self._pi[1], self._pi[2])

#---------------------------------------------------------------------------------------------------------------------- 
cdef class PairedSnvMixParameters(Parameters):
    def __init__(self, **kwargs):
        default_mu = (
                      (100, 2),
                      (50, 50),
                      (2, 100)
                      )
        
        default_pi = (1000, 100, 100)
        
        mu_N = kwargs.get('mu_N', default_mu)
        mu_T = kwargs.get('mu_T', default_mu)

        pi_N = kwargs.get('pi_N', default_pi)
        pi_T = kwargs.get('pi_T', default_pi)
        
        self._normal_priors = SnvMixPriors(mu=mu_N, pi=pi_N)
        self._tumour_priors = SnvMixPriors(mu=mu_T, pi=pi_T)
    
    property normal:
        def __get__(self):
            return self._normal_priors
    
    property tumour:
        def __get__(self):
            return self._tumour_priors

#=======================================================================================================================
# Parameters
#=======================================================================================================================
cdef class SnvMixParameters(Parameters):
    def __init__(self, **kwargs):        
        self._priors = kwargs.get('priors', SnvMixPriors())
        
        mu = kwargs.get('mu', (0.99, 0.5, 0.01))
        pi = kwargs.get('pi', (0.99, 0.009, 0.001))
        
        for g in range(NUM_GENOTYPES):
            self._mu[g] = mu[g]
            self._pi[g] = pi[g]
    
    property mu:
        def __get__(self):
            return tuple([x for x in self._mu[:NUM_GENOTYPES]])
    
    property pi:
        def __get__(self):
            return tuple([x for x in self._pi[:NUM_GENOTYPES]])

    cdef _normalise_pi(self):
        cdef int g
        cdef double norm_const
        
        norm_const = 0
        
        for g in range(NUM_GENOTYPES):
            norm_const += self._pi[g]
        
        for g in range(NUM_GENOTYPES):
            self._pi[g] = self._pi[g] / norm_const

    cdef update(self, double * n, double * a, double * b):
        self._update_mu(a, b)
        self._update_pi(n)

    cdef _update_mu(self, double * a, double * b):
        cdef int g
        cdef double alpha, beta, denom
    
        for g in range(NUM_GENOTYPES):
            alpha = a[g] + self._priors._mu[g][0] - 1
            beta = b[g] + self._priors._mu[g][1] - 1
            denom = alpha + beta

            self._mu[g] = alpha / denom
            
    cdef _update_pi(self, double * n):
        cdef int g
        cdef double denom
        cdef double pi[NUM_GENOTYPES]
        
        denom = 0
        
        for g in range(NUM_GENOTYPES):
            pi[g] = n[g] + self._priors._pi[g] - 1
            denom += pi[g]
        
        for g in range(NUM_GENOTYPES):
            self._pi[g] = pi[g] / denom
    
    cdef double _get_prior_log_likelihood(self):
        cdef double ll
        cdef double x[2]
        
        ll = 0
        
        for g in range(NUM_GENOTYPES):
            x[0] = self._mu[g]
            x[1] = 1 - self._mu[g]
            
            ll += dirichlet_log_likelihood(x, self._priors._mu[g], NUM_BASES)
        
        ll += dirichlet_log_likelihood(self._pi, self._priors._pi, NUM_GENOTYPES)
        
        return ll
            
#---------------------------------------------------------------------------------------------------------------------- 
cdef class PairedSnvMixParameters(object):
    def __init__(self, **kwargs):
        cdef PairedSnvMixPriors priors
        cdef SnvMixPriors normal_priors, tumour_priors
        
        priors = kwargs.get('priors', PairedSnvMixPriors())
        normal_priors = priors._normal_priors
        tumour_priors = priors._tumour_priors
         
        mu_N = kwargs.get('mu_N', (0.99, 0.5, 0.01))
        mu_T = kwargs.get('mu_T', (0.99, 0.5, 0.01))
        pi_N = kwargs.get('pi_N', (0.99, 0.009, 0.001))
        pi_T = kwargs.get('pi_T', (0.99, 0.009, 0.001))
        
        self._normal_params = SnvMixParameters(priors=normal_priors, mu=mu_N, pi=pi_N)
        self._tumour_params = SnvMixParameters(priors=tumour_priors, mu=mu_T, pi=pi_T)
    
    property normal:
        def __get__(self):
            return self._normal_params
    
    property tumour:
        def __get__(self):
            return self._tumour_params
    
    property mu_N:
        def __get__(self):
            return self._normal_params.mu

    property mu_T:
        def __get__(self):
            return self._tumour_params.mu
        
    property pi_N:
        def __get__(self):
            return self._normal_params.pi

    property pi_T:
        def __get__(self):
            return self._tumour_params.pi 
