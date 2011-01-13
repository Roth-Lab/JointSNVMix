'''
Created on 2010-11-17

@author: andrew
'''
import numpy as np

from scipy.special import psi

class _Posterior( object ):
    def __init__( self, priors ):
        self.priors = priors
        self.parameters = _Parameters()
                
    def set_latent_variables( self, latent_variables ):
        self.latent_variables = latent_variables

    def update( self ):
        N_g = self.latent_variables.responsibilities.sum( axis=0 )
        
        self.kappa_bar = self.priors.kappa + N_g
        
        self.expected_log_pi = psi( self.kappa_bar ) - psi( self.kappa_bar.sum() )
        
        self.alpha_bar = []
        self.beta_bar = []
        self.expected_log_mu = []
        self.expected_log_one_minus_mu = []
        
        for genome in range( 2 ):
            self.alpha_bar.append( self.latent_variables.marginals[genome] * \
                                   self.latent_variables.a_plus_alpha[genome] + \
                                   1 )
            
            self.beta_bar.append( self.latent_variables.marginals[genome] * \
                                  self.latent_variables.b_plus_beta[genome] + \
                                  1 )
            
            self.expected_log_mu.append( psi( self.alpha_bar[genome] ) - \
                                         psi( self.alpha_bar[genome] + self.beta_bar[genome] ) )
            
            self.expected_log_one_minus_mu.append( psi( self.beta_bar[genome] ) - \
                                                   psi( self.alpha_bar[genome] + self.beta_bar[genome] ) )

class _Parameters( object ):
    def __init__( self ):
        self.alpha = []
        self.beta = []
        
        for i in range( 2 ):
            self.alpha.append( np.array( [99., 5., 1.] ) )
            self.beta.append( np.array( [1., 5., 99.] ) )

