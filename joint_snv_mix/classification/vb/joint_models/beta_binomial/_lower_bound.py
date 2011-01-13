'''
Created on 2010-11-17

@author: Andrew Roth
'''
import numpy as np

from scipy.special import betaln
from ....utils.log_pdf import log_dirichlet_constant

class _LowerBound( object ):
    def __init__( self, priors, latent_variables, posterior ):
        self.priors = priors
        self.lv = latent_variables
        self.posterior = posterior
        
    def get_lower_bound( self ):
        lower_bound = 0
        
#        lower_bound += self._get_expectation_log_p_x_plus_log_mu()
        lower_bound += self._get_expectation_log_p_z()
        lower_bound += self._get_expectation_log_p_pi()
        
        lower_bound -= self._get_expectation_log_q_z()
        lower_bound -= self._get_expectation_log_q_pi()
#        lower_bound -= self._get_expectation_log_q_mu()

        for i in range( 2 ):
            lower_bound += np.sum( betaln( self.posterior.alpha_bar[i], self.posterior.beta_bar[i] ) - \
                                   betaln( self.posterior.parameters.alpha[i], self.posterior.parameters.beta[i] ) )
        
        return lower_bound

    def _get_expectation_log_p_x_plus_log_mu( self ):
        terms = []
        
        for i in range( 2 ):
            term = ( self.posterior.alpha_bar[i] - 1 ) * self.posterior.expected_log_mu[i] + \
                   ( self.posterior.beta_bar[i] - 1 ) * self.posterior.expected_log_one_minus_mu[i] 
            term = term.sum()
            
            terms.append( term )
        
        return sum( terms )
    
    def _get_expectation_log_p_z( self ):
        term = self.posterior.expected_log_pi * self.lv.responsibilities.sum( axis=0 )
        term = term.sum()
        
        return term
    
    def _get_expectation_log_p_pi( self ):
        term_1 = log_dirichlet_constant( self.priors.kappa )
        
        term_2 = ( self.priors.kappa - 1. ) * self.posterior.expected_log_pi
        term_2 = term_2.sum()
        
        return term_1 + term_2

    
    def _get_expectation_log_q_z( self ):
        eps = np.finfo( np.float64 ).eps
        
        r = self.lv.responsibilities + eps
        
        term = r * np.log( r )
        term = term.sum( axis=1 )
        term = term.sum()
        
        return term
    
    def _get_expectation_log_q_pi( self ):
        k = self.posterior.kappa_bar
        
        term_1 = log_dirichlet_constant( k )
        
        term_2 = ( k - 1 ) * self.posterior.expected_log_pi
        term_2 = term_2.sum()
        
        return term_1 + term_2
    
    def _get_expectation_log_q_mu( self ):
        alpha_bar = self.posterior.alpha_bar
        beta_bar = self.posterior.beta_bar
        
        terms = []
        
        for i in range( 2 ):
            term = betaln( alpha_bar[i] , beta_bar[i] ) + \
                   ( alpha_bar[i] - 1 ) * self.posterior.expected_log_mu[i] + \
                   ( beta_bar[i] - 1 ) * self.posterior.expected_log_one_minus_mu[i]
            
            term = term.sum()
            
            terms.append( term )
        
        return sum( terms )
        
        
