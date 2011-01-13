'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from jsm_models.utils.log_pdf import log_dirichlet_pdf
 
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
