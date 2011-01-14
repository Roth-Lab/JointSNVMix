'''
Created on 2010-11-17

@author: andrew
'''
from scipy.special import psi

class _Posterior( object ):
    def __init__( self, Priors ):
        self.priors = Priors
        
        self.kappa_bar = self.priors['kappa']
        
        self.alpha_bar_1 = self.priors['alpha'][0]
        self.alpha_bar_2 = self.priors['alpha'][1]
        
        self.beta_bar_1 = self.priors['beta'][0]
        self.beta_bar_2 = self.priors['beta'][1]
        
        self._update_auxillary_statistics()
    
    def update( self, ESS ):
        self.kappa_bar = self.priors['kappa'] + ESS.N_g
        
        self.alpha_bar_1 = self.priors['alpha'][0] + ESS.a_bar_1        
        self.alpha_bar_2 = self.priors['alpha'][1] + ESS.a_bar_2
        
        self.beta_bar_1 = self.priors['beta'][0] + ESS.b_bar_1     
        self.beta_bar_2 = self.priors['beta'][1] + ESS.b_bar_2
        
        self._update_auxillary_statistics()

    def _update_auxillary_statistics( self ):
        self.expected_log_pi = psi( self.kappa_bar ) - psi( self.kappa_bar.sum() )
        
        self.expected_log_mu_1 = psi( self.alpha_bar_1 ) - psi( self.alpha_bar_1 + self.beta_bar_1 )
        self.expected_log_mu_2 = psi( self.alpha_bar_2 ) - psi( self.alpha_bar_2 + self.beta_bar_2 )

        self.expected_log_one_minus_mu_1 = psi( self.beta_bar_1 ) - psi( self.alpha_bar_1 + self.beta_bar_1 )
        self.expected_log_one_minus_mu_2 = psi( self.beta_bar_2 ) - psi( self.alpha_bar_2 + self.beta_bar_2 )
