from joint_snv_mix.classification.base import *
from joint_snv_mix.classification.joint import *

from joint_snv_mix.classification.utils.log_pdf import log_beta_pdf, log_binomial_likelihood, get_joint_log_likelihoods

class JointBinomialRunner(ModelRunner):
    def __init__(self):
        self.model = JointBinomialModel()
        self.priors_parser = JointBinomialPriorParser()
        self.parameter_parser = JointBinomialParameterParser()

class JointBinomialModel(EMModel):
    def __init__(self):
        self.trainer_class = JointBinomialModelTrainer
        
        self.log_likelihood_func = joint_binomial_log_likelihood

class JointBinomialModelTrainer(EMModelTrainer):
    def _init_components(self):
        self.latent_variables = JointBinomialLatentVariables(self.data)
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = JointBinomialPosterior(self.data, self.priors, self.responsibilities)
        
        self.lower_bound = JointBinomialLowerBound(self.data, self.priors)

class JointBinomialLatentVariables(JointLatentVariables):
    def __init__(self, data):
        JointLatentVariables.__init__(self, data)
        
        self.likelihood_func = joint_binomial_log_likelihood

class JointBinomialPosterior(EMPosterior):
    def __init__(self, data, priors, responsibilities, nclass=3):
        self.nclass = nclass
        
        EMPosterior.__init__(self, data, priors, responsibilities)
    
    def _init_parameters(self):
        self.parameters = {}
        self.parameters['normal'] = {}
        self.parameters['tumour'] = {}
        
        self._update_mix_weights()
        
        self._update_density_parameters()
           
    def _update_density_parameters(self):
        marginals = get_marginals(self.responsibilities, self.nclass)
        
        for genome in constants.genomes:
            a = self.data.a[genome]
            b = self.data.b[genome]
            
            alpha = self.priors[genome]['mu']['alpha']
            beta = self.priors[genome]['mu']['beta']
            
            tau = marginals[genome]
            
            self.parameters[genome]['mu'] = self._update_mu(a, b, alpha, beta, tau)
    
    def _update_mu(self, a, b, alpha, beta, tau):       
        d = a + b
        
        n = self.data.nrows
        shape = (n, 1)
        
        a = a.reshape(shape)
        d = d.reshape(shape)
        
        ref_sum = np.sum(tau * a, axis=0)
        
        depth_sum = np.sum(tau * d, axis=0)
        
        numerator = ref_sum + alpha - 1.
        
        denominator = depth_sum + alpha + beta - 2.
        
        return np.exp(np.log(numerator) - np.log(denominator))

class JointBinomialLowerBound(EMLowerBound):
    def __init__(self, data, priors):  
        EMLowerBound.__init__(self, data, priors)
        
        self.log_likelihood_func = joint_binomial_log_likelihood
    
    def _get_log_density_parameters_prior(self):
        log_prior = 0.
        
        for genome in constants.genomes:
            mu = self.parameters[genome]['mu']
            
            alpha = self.priors[genome]['mu']['alpha']
            beta = self.priors[genome]['mu']['beta']
            
            log_prior += log_beta_pdf(mu, alpha, beta)
        
        log_prior = log_prior.sum()

        return log_prior
    
class JointBinomialPriorParser(JointModelPriorParser):
    def __init__(self):
        JointModelPriorParser.__init__(self)
        
        self.parameter_names = ('mu',)
        
        self.hyper_parameter_names = {}
        self.hyper_parameter_names['mu'] = ('alpha', 'beta')

class JointBinomialParameterParser(JointParameterParser):
    def __init__(self):
        JointParameterParser.__init__(self)
        
        self.parameter_names = ('mu',)

def get_marginals(responsibilities, nclass):
    nrows = responsibilities.shape[0]

    shape = (nrows, nclass, nclass)
    
    responsibilities = responsibilities.reshape(shape)
    
    marginals = {}
    
    marginals['normal'] = responsibilities.sum(axis=2)
    marginals['tumour'] = responsibilities.sum(axis=1)

    return marginals

def joint_binomial_log_likelihood(data, parameters):
    log_likelihoods = {}
    
    for genome in constants.genomes:
        a = data.a[genome]
        b = data.b[genome]
        d = a + b
        
        mu = parameters[genome]['mu']
    
        log_likelihoods[genome] = log_binomial_likelihood(a, d, mu)

    pi = parameters['pi']

    log_likelihoods = get_joint_log_likelihoods(log_likelihoods, pi)
    
    return log_likelihoods
