'''
Created on 2010-08-30

@author: Andrew
'''
import numpy as np

from scipy.special import betaln, gammaln

def log_beta_pdf( x, a, b ):
    log_p = -betaln( a, b ) + ( a - 1 ) * np.log( x ) + ( b - 1 ) * np.log( 1 - x )

    return log_p

def log_dirichlet_pdf( x, kappa ):
    if  np.abs( x.sum() - 1. ) > 1e-12:
        return np.float( '-inf' )
    
    log_normalisation_constant = gammaln( kappa.sum() ) - gammaln( kappa ).sum()

    log_likelihood = np.sum( ( kappa - 1 ) * np.log( x ) )

    log_p = log_normalisation_constant + log_likelihood

    return log_p

def log_gamma_pdf( x, shape, scale ):   
    log_p = ( shape - 1 ) * np.log( x ) - gammaln( shape ) - shape * np.log( scale ) - ( x / scale )
    
    return log_p

def log_beta_binomial_likelihood( k, n, alpha, beta ):
    column_shape = ( k.size, 1 )
    k = k.reshape( column_shape )
    n = n.reshape( column_shape )

    row_shape = ( 1, alpha.size )
    alpha = alpha.reshape( row_shape )
    beta = beta.reshape( row_shape )

    return betaln( k + alpha, n - k + beta ) - betaln( alpha, beta )

def log_dirichlet_constant( kappa ):
    return gammaln( kappa.sum() ) - gammaln( kappa ).sum()

def log_binomial_likelihood( k, n, mu ):
    column_shape = ( k.size, 1 )
    k = k.reshape( column_shape )
    n = n.reshape( column_shape )
    
    row_shape = ( 1, mu.size )
    mu = mu.reshape( row_shape )
    
    return k * np.log( mu ) + ( n - k ) * np.log( 1 - mu )
