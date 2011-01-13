'''
Created on 2010-11-17

@author: Andrew Roth
'''
import numpy as np

from scipy.special import betaln
from ....utils.log_pdf import log_dirichlet_constant

class _LowerBound( object ):
    def __init__( self, Priors, Posterior, ESS, LatentVariables ):
        self.priors = Priors
        self.posterior = Posterior
        self.ess = ESS
        self.lv = LatentVariables
        
    def get_lower_bound( self ):
        lower_bound = 0
        
        lower_bound += self._get_expectation_log_p_x()
        lower_bound += self._get_expectation_log_p_z()
        lower_bound += self._get_expectation_log_p_pi()
        lower_bound += self._get_expectation_log_p_mu()
        
        lower_bound -= self._get_expectation_log_q_z()
        lower_bound -= self._get_expectation_log_q_pi()
        lower_bound -= self._get_expectation_log_q_mu()
        
        return lower_bound

    def _get_expectation_log_p_x( self ):
        term_1 = self.ess.a_bar_1 * self.posterior.expected_log_mu_1 + \
                 self.ess.b_bar_1 * self.posterior.expected_log_one_minus_mu_1 
        term_1 = term_1.sum()
        
        term_2 = self.ess.a_bar_2 * self.posterior.expected_log_mu_2 + \
                 self.ess.b_bar_2 * self.posterior.expected_log_one_minus_mu_2 
        term_2 = term_2.sum()
        
        return term_1 + term_2
    
    def _get_expectation_log_p_z( self ):
        term = self.posterior.expected_log_pi * self.ess.N_g
        term = term.sum()
        
        return term
    
    def _get_expectation_log_p_pi( self ):
        term_1 = log_dirichlet_constant( self.priors.kappa )
        
        term_2 = ( self.priors.kappa - 1. ) * self.posterior.expected_log_pi
        term_2 = term_2.sum()
        
        return term_1 + term_2

    def _get_expectation_log_p_mu( self ):
        alpha_1 = self.priors.alpha_1
        beta_1 = self.priors.beta_1
        
        term_1 = -betaln( alpha_1, beta_1 ) + \
                 ( alpha_1 - 1 ) * self.posterior.expected_log_mu_1 + \
                 ( beta_1 - 1 ) * self.posterior.expected_log_one_minus_mu_1
        term_1 = term_1.sum()
        
        alpha_2 = self.priors.alpha_2
        beta_2 = self.priors.beta_2
        
        term_2 = -betaln( alpha_2, beta_2 ) + \
                ( alpha_2 - 1 ) * self.posterior.expected_log_mu_2 + \
                ( beta_2 - 1 ) * self.posterior.expected_log_one_minus_mu_2
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
        alpha_1 = self.posterior.alpha_bar_1
        beta_1 = self.posterior.beta_bar_1
        
        term_1 = -betaln( alpha_1, beta_1 ) + \
                 ( alpha_1 - 1 ) * self.posterior.expected_log_mu_1 + \
                 ( beta_1 - 1 ) * self.posterior.expected_log_one_minus_mu_1
        term_1 = term_1.sum()
        
        alpha_2 = self.posterior.alpha_bar_2
        beta_2 = self.posterior.beta_bar_2
        
        term_2 = -betaln( alpha_2, beta_2 ) + \
                ( alpha_2 - 1 ) * self.posterior.expected_log_mu_2 + \
                ( beta_2 - 1 ) * self.posterior.expected_log_one_minus_mu_2
        term_2 = term_2.sum()
        
        return term_1 + term_2
        
        
