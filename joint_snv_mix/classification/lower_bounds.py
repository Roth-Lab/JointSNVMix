'''
Created on 2011-01-18

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix.classification.utils.log_pdf import log_dirichlet_pdf, log_beta_pdf, log_gamma_pdf, \
    log_translated_gamma_pdf
from joint_snv_mix.classification.likelihoods import independent_binomial_log_likelihood, \
    independent_beta_binomial_log_likelihood, joint_beta_binomial_log_likelihood, joint_binomial_log_likelihood
from joint_snv_mix import constants

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
        alpha = self.parameters['alpha']
        beta = self.parameters['beta']
            
        s = alpha + beta
        mu = alpha / s
    
        precision_priors = self.priors['precision']
        location_priors = self.priors['location']
            
        precision_term = np.sum( log_gamma_pdf( s, precision_priors['shape'], precision_priors['scale'] ) )
        location_term = np.sum( log_beta_pdf( mu, location_priors['alpha'], location_priors['beta'] ) )
                
        return precision_term + location_term
    
class IndependenBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = independent_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        log_prior = 0.
            
        mu = self.parameters['mu']
        
        alpha = self.priors['mu']['alpha']
        beta = self.priors['mu']['beta']
        
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
        precision_term = 0
        location_term = 0
        
        genomes = ['junk']
        genomes.extend( constants.genomes )
        
        for genome in genomes:   
            alpha = self.parameters[genome]['alpha']
            beta = self.parameters[genome]['beta']
            
            s = alpha + beta
            mu = alpha / s
    
            precision_priors = self.priors[genome]['precision']
            location_priors = self.priors[genome]['location']
            
            precision_term += np.sum( log_translated_gamma_pdf( s,
                                                                precision_priors['shape'],
                                                                precision_priors['scale'],
                                                                precision_priors['min'] 
                                                                ) )
            
            location_term += np.sum( log_beta_pdf( mu,
                                                   location_priors['alpha'],
                                                   location_priors['beta'] ) )
                    
        return precision_term + location_term

class JointBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        log_prior = 0.
        
        for genome in constants.genomes:
            mu = self.parameters[genome]['mu']
            
            alpha = self.priors[genome]['mu']['alpha']
            beta = self.priors[genome]['mu']['beta']
            
            log_prior += log_beta_pdf( mu, alpha, beta )
        
        log_prior = log_prior.sum()

        return log_prior
