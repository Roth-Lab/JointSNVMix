'''
The independent binomial paired data classifier. This is equivalent to using the SNVMix1 model on the tumour normal pair
and then multiplying the probabilities to get the joint genotype probabilities.


Created on 2011-03-31

@author: Andrew Roth
'''
class IndependentBinomialRunner(IndependentModelRunner):
    def __init__(self):
        self.model = IndependentBinomialModel()
        self.priors_parser = IndependentBinomialPriorParser()
        self.parameter_parser = IndependentBinomialParameterParser()

class IndependentBinomialModel(EMModel):
    def __init__(self):
        self.trainer_class = IndependentBinomialModelTrainer
        
        self.log_likelihood_func = independent_binomial_log_likelihood

class IndependentBinomialModelTrainer(EMModelTrainer):
    def _init_components(self):
        self.latent_variables = IndependentBinomialLatentVariables(self.data)
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = IndependentBinomialPosterior(self.data, self.priors, self.responsibilities)
        
        self.lower_bound = IndependenBinomialLowerBound(self.data, self.priors)

class IndependentBinomialLatentVariables(IndependentLatenVariables):
    def __init__(self, data):
        IndependentLatenVariables.__init__(self, data)
        
        self.likelihood_func = independent_binomial_log_likelihood
        
class IndependentBinomialPosterior(EMPosterior):
    def _init_parameters(self):
        self.parameters = {}
        
        self._update_mix_weights()

        alpha = self.priors['mu']['alpha']
        beta = self.priors['mu']['beta']
        
        self.parameters['mu'] = alpha / (alpha + beta)
        
        print self.parameters

    def _update_density_parameters(self):
        a = self.data.a
        b = self.data.b
        d = a + b
        
        shape = (self.data.nrows, 1)
        a = a.reshape(shape)
        d = d.reshape(shape)
        
        alpha = self.priors['mu']['alpha']
        beta = self.priors['mu']['beta']
        
        resp = self.responsibilities
        
        a_bar = np.sum(resp * a, axis=0)
        
        d_bar = np.sum(resp * d, axis=0)
        
        numerator = a_bar + alpha - 1
        denominator = d_bar + alpha + beta - 2
        
        self.parameters['mu'] = numerator / denominator
        
class IndependenBinomialLowerBound(EMLowerBound):
    def __init__(self, data, priors):  
        EMLowerBound.__init__(self, data, priors)
        
        self.log_likelihood_func = independent_binomial_log_likelihood
    
    def _get_log_density_parameters_prior(self):
        log_prior = 0.
            
        mu = self.parameters['mu']
        
        alpha = self.priors['mu']['alpha']
        beta = self.priors['mu']['beta']
        
        log_prior += log_beta_pdf(mu, alpha, beta)
    
        log_prior = log_prior.sum()

        return log_prior
           
class IndependentBinomialPriorParser(IndepedendentPriorParser):
    def __init__(self):
        IndepedendentPriorParser.__init__(self)
        
        self.parameter_names = ('mu',)
        
        self.hyper_parameter_names = {}
        self.hyper_parameter_names['mu'] = ('alpha', 'beta')

class IndependentBinomialParameterParser(IndepedendentParameterParser):
    def __init__(self):
        IndepedendentParameterParser.__init__(self)
        
        self.parameter_names = ('mu',)
        
def independent_binomial_log_likelihood(data, parameters):
    a = data.a
    b = data.b
    
    d = a + b
    
    mu = parameters['mu']

    log_likelihoods = log_binomial_likelihood(a, d, mu)

    pi = parameters['pi']
    log_pi = np.log(pi)

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods

