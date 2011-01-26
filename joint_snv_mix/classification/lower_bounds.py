'''
Created on 2011-01-18

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix.classification.utils.log_pdf import log_dirichlet_pdf, log_beta_pdf, log_gamma_pdf, \
    log_translated_gamma_pdf
from joint_snv_mix.classification.likelihoods import independent_binomial_log_likelihood, \
    independent_beta_binomial_log_likelihood, joint_beta_binomial_log_likelihood, joint_binomial_log_likelihood

#=======================================================================================================================
# Abstract Models
#=======================================================================================================================
class EMLowerBound( object ):
    def __init__( self, data, priors ):                
        self.priors = priors
        self.data = data
        
        self.log_likelihood_func = None

    def get_lower_bound( self, parameters ):
        self.parameters = parameters
        
        log_likelihood = self._get_log_likelihood()

        log_mix_weight_prior = self._get_log_mix_weight_prior()

        log_density_parameters_prior = self._get_log_density_parameters_prior()

        lower_bound = log_likelihood + log_mix_weight_prior + log_density_parameters_prior

        return lower_bound
    
    def _get_log_likelihood( self ):
        log_likelihoods = self.log_likelihood_func( self.data, self.parameters )
        
        log_likelihood = np.logaddexp.reduce( log_likelihoods, axis=1 ).sum()

        return log_likelihood

    def _get_log_mix_weight_prior( self ):
        pi = self.parameters['pi']
        kappa = self.priors['kappa']
        
        log_prior = log_dirichlet_pdf( pi, kappa )
        log_prior = log_prior.sum()

        return log_prior

    def _get_log_density_parameters_prior( self ):
        raise NotImplemented
    
#=======================================================================================================================
# Independent Models
#=======================================================================================================================
class IndependentBetaBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = independent_beta_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        precision_term = 0
        location_term = 0
        
        for component in range( 3 ):
            alpha = self.parameters['alpha'][component]
            beta = self.parameters['beta'][component]
            
            s = alpha + beta
            mu = alpha / s
    
            precision_priors = self.priors['precision'][component]
            location_priors = self.priors['location'][component]
            
            precision_term += log_gamma_pdf( s, precision_priors[0], precision_priors[1] )
            location_term += log_beta_pdf( mu, location_priors[0], location_priors[1] )
                
        return precision_term + location_term
    
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
    
#=======================================================================================================================
# Joint Models
#=======================================================================================================================
class JointBetaBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_beta_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        alpha_term = 0
        beta_term = 0
        
        for genome in range( 2 ):
            for component in range( 3 ):
                alpha = self.parameters['alpha'][genome][component]
                beta = self.parameters['beta'][genome][component]
                
                alpha_priors = self.priors['alpha'][genome][component]
                beta_priors = self.priors['beta'][genome][component]
                
                alpha_term += log_translated_gamma_pdf( alpha, alpha_priors[0], alpha_priors[1], alpha_priors[2] )
                beta_term += log_translated_gamma_pdf( beta, beta_priors[0], beta_priors[1], alpha_priors[2] )
                
        return alpha_term + beta_term

class JointBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        log_prior = 0.
        
        for genome in range( 2 ):
            mu = self.parameters['mu'][genome]
            alpha = self.priors['alpha'][genome]
            beta = self.priors['beta'][genome]
            
            log_prior += log_beta_pdf( mu, alpha, beta )
        
        log_prior = log_prior.sum()

        return log_prior
