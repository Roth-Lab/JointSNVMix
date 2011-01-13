'''
Created on 2010-11-17

@author: Andrew Roth
'''
import numpy as np

from scipy.special import betaln

from ....utils.marginal_responsibilities import get_marginals
from ....utils.normalise import log_space_normalise_rows 

class _LatentVariables( object ):
    def __init__( self, data, priors ):
        self.data = data
        self.priors = priors
        
    def set_posterior( self, posterior ):
        self.posterior = posterior
        self._add_priors_to_data()
        
    def set_responsibilities( self, repsonsibilities ):
        self.responsibilities = repsonsibilities
        self.marginals = get_marginals( self.responsibilities )

    def update( self ):       
        alpha = self.posterior.parameters.alpha
        beta = self.posterior.parameters.beta
        
        log_mu = self.posterior.expected_log_mu
        log_one_minus_mu = self.posterior.expected_log_one_minus_mu
        
        log_pi = self.posterior.expected_log_pi
        
        shape = ( 2, self.data.nrows, 3 )      
        term = np.zeros( shape )
        
        for i in range( 2 ):            
            term[i, :, :] = self.a_plus_alpha[i] * log_mu[i] + \
                            self.b_plus_beta[i] * log_one_minus_mu[i] - \
                            betaln( alpha[i], beta[i] )

        columns = []

        for i in range( 3 ):
            for j in range( 3 ):
                k = 3 * i + j
                
                columns.append( log_pi[k] + term[0, :, i] + term[1, :, j] )
                
        self.responsibilities = np.column_stack( columns )

        self.responsibilities = log_space_normalise_rows( self.responsibilities )
        
        self.marginals = get_marginals( self.responsibilities )
        
    def _add_priors_to_data( self ):
        n = self.data.nrows
        shape = ( n, 1 )
        
        self.a_plus_alpha = []
        self.b_plus_beta = []
        
        for i in range( 2 ):
            a = self.data.a[i].reshape( shape )
            b = self.data.b[i].reshape( shape )
            
            alpha = self.posterior.parameters.alpha[i]
            beta = self.posterior.parameters.beta[i]
            
            self.a_plus_alpha.append( a + alpha - 1 )
            self.b_plus_beta.append( b + beta - 1 )




