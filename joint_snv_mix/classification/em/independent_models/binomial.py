'''
Created on 2010-12-09

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix.classification.em.em_lower_bound import EMLowerBound
from joint_snv_mix.classification.em.em_model import EMModel, EMModelTrainer
from joint_snv_mix.classification.em.em_posterior import EMPosterior
from joint_snv_mix.classification.em.independent_models.independent_latent_variables import IndependentLatenVariables
from joint_snv_mix.classification.utils.log_pdf import log_beta_pdf, log_binomial_likelihood

class IndependentBinomialModel( EMModel ):
    def __init__( self ):
        self.trainer_class = IndependentBinomialModelTrainer
        
        self.log_likelihood_func = independent_binomial_log_likelihood

class IndependentBinomialModelTrainer( EMModelTrainer ):
    def _init_components( self ):
        self.latent_variables = IndependentBinomialLatentVariables( self.data )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = IndependentBinomialPosterior( self.data, self.priors, self.responsibilities )
        
        self.lower_bound = IndependenBinomialLowerBound( self.data, self.priors )
        
class IndependentBinomialLatentVariables( IndependentLatenVariables ):
    def __init__( self, data ):
        IndependentLatenVariables.__init__( self, data )
        
        self.likelihood_func = independent_binomial_log_likelihood
        
class IndependenBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = independent_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        log_prior = 0.
            
        mu = self.parameters['mu']
        alpha = self.priors['alpha']
        beta = self.priors['beta']
        
        log_prior += log_beta_pdf( mu, alpha, beta )
    
        log_prior = log_prior.sum()

        return log_prior
    
class IndependentBinomialPosterior( EMPosterior ):
    def _init_parameters( self ):
        self.parameters = {}
        
        self._update_mix_weights()

        alpha = self.priors['alpha']
        beta = self.priors['beta']
        
        self.parameters['mu'] = alpha / ( alpha + beta )
        
        print self.parameters

    def _update_density_parameters( self ):
        a = self.data.a
        b = self.data.b
        d = a + b
        
        shape = ( self.data.nrows, 1 )
        a = a.reshape( shape )
        d = d.reshape( shape )
        
        alpha = self.priors['alpha']
        beta = self.priors['beta']
        
        resp = self.responsibilities
        
        a_bar = np.sum( resp * a, axis=0 )
        
        d_bar = np.sum( resp * d, axis=0 )
        
        numerator = a_bar + alpha - 1
        denominator = d_bar + alpha + beta - 2
        
        self.parameters['mu'] = numerator / denominator

def independent_binomial_log_likelihood( data, parameters ):
    a = data.a
    b = data.b
    
    d = a + b
    
    mu = parameters['mu']

    log_likelihoods = log_binomial_likelihood( a, d, mu )

    pi = parameters['pi']
    log_pi = np.log( pi )

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods
