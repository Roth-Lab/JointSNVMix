'''
Created on 2010-12-09

@author: Andrew Roth
'''
import multiprocessing

import numpy as np

from joint_snv_mix.classification.em.em_lower_bound import EMLowerBound
from joint_snv_mix.classification.em.em_model import EMModel, EMModelTrainer
from joint_snv_mix.classification.em.em_posterior import EMPosterior
from joint_snv_mix.classification.em.independent_models.independent_latent_variables import IndependentLatenVariables
from joint_snv_mix.classification.utils.beta_binomial_map_estimators import get_mle_p
from joint_snv_mix.classification.utils.log_pdf import log_gamma_pdf, log_beta_pdf, log_beta_binomial_likelihood



class IndependenBetaBinomialModel( EMModel ):
    def __init__( self ):
        self.trainer_class = IndependenBetaBinomialTrainer
        
        self.log_likelihood_func = independent_beta_binomial_log_likelihood

class IndependenBetaBinomialTrainer( EMModelTrainer ):
    def _init_components( self ):
        self.latent_variables = IndependentBetaBinomialLatentVariables( self.data )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = IndependentBetaBinomialPosterior( self.data, self.priors, self.responsibilities )
        
        self.lower_bound = IndependentBetaBinomialLowerBound( self.data, self.priors )

class IndependentBetaBinomialLatentVariables( IndependentLatenVariables ):
    def __init__( self, data ):
        IndependentLatenVariables.__init__( self, data )
        
        self.likelihood_func = independent_beta_binomial_log_likelihood
        
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
    
class IndependentBetaBinomialPosterior( EMPosterior ):
    def __init__( self, *args ):
        EMPosterior.__init__( self, *args )
        
        self.p = multiprocessing.Pool( processes=3, maxtasksperchild=1 )
    
    def _init_parameters( self ):
        '''
        Initialise parameters using method of moments (MOM) estiamtes.
        '''
        self.parameters = {}
        
        self._update_mix_weights()
        
#        s = self.priors['precision'][:, 0] * self.priors['precision'][ :, 1] 
#        
#        mu = self.priors['location'][:, 0] / \
#            ( self.priors['location'][:, 0] + self.priors['location'][ :, 1] ) 
#        
#        self.parameters['alpha'] = s * mu
#        self.parameters['beta'] = s * ( 1 - mu )

        self.parameters['alpha'] = np.array( [99, 5, 1], dtype=np.float64 )
        self.parameters['beta'] = np.array( [1, 5, 99], dtype=np.float64 )
        
        print "Initial parameter values : ", self.parameters
    
    def _update_density_parameters( self ):        
        a = self.data.a
        b = self.data.b
        
        vars = []
              
        for component in range( 3 ):
            x = np.zeros( ( 2, ) )
            
            x[0] = self.parameters['alpha'][component]
            x[1] = self.parameters['beta'][component]
            
            resp = self.responsibilities[:, component]
            
            precision_prior = self.priors['precision'][component]
            location_prior = self.priors['location'][component]
            
            vars.append( [x, a, b, resp, component, precision_prior, location_prior] )
                
        
        results = self.p.map( get_mle_p, vars )
        
        for component in range( 3 ):
            self.parameters['alpha'][component] = results[component][0]
            self.parameters['beta'][component] = results[component][1]
            
def independent_beta_binomial_log_likelihood( data, parameters ):
    a = data.a
    b = data.b
    
    d = a + b
    
    alpha = parameters['alpha']
    beta = parameters['beta']
    
    log_likelihoods = log_beta_binomial_likelihood( a, d, alpha, beta )

    pi = parameters['pi']
    log_likelihoods = log_likelihoods + np.log( pi )
    
    return log_likelihoods
