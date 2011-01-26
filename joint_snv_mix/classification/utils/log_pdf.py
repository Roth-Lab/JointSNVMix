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
#    if  np.abs( x.sum() - 1. ) > 1e-12:
#        return np.float( '-inf' )
    if kappa.ndim == 1: 
        log_normalisation_constant = gammaln( kappa.sum() ) - gammaln( kappa ).sum()
    
        log_likelihood = np.sum( ( kappa - 1 ) * np.log( x ) )
    else:
        log_normalisation_constant = gammaln( kappa.sum( axis=1 ) ) - gammaln( kappa ).sum( axis=1 )
    
        log_likelihood = np.sum( ( kappa - 1 ) * np.log( x ), axis=1 )

    log_p = log_normalisation_constant + log_likelihood

    return log_p

def log_gamma_pdf( x, shape, scale ):
    print x, shape, scale
    
    log_p = ( shape - 1 ) * np.log( x ) - gammaln( shape ) - shape * np.log( scale ) - ( x / scale )
    
    return log_p

def log_translated_gamma_pdf( x, scale, shape, min ):
    x = ( x - min )
    
    return log_gamma_pdf( x, shape, scale )

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

def log_multinomial_likelihood( counts, rho ):
    '''
    counts is Nxk and rho is cxk
    
    p is Nxc
    '''
    
    counts = counts.reshape( ( counts.shape[0], counts.shape[1], 1 ) )
    
    rho = np.swapaxes( rho, 0, 1 )
    rho = rho.reshape( ( 1, rho.shape[0], rho.shape[1] ) )
    
    p = counts * np.log( rho )
    
    p = p.sum( axis=1 )
    
    return p

def log_multivariate_polya_likelihood( counts, alpha ):
    '''
    counts is Nxk and rho is cxk
    
    p is Nxc
    '''
    counts = counts.reshape( ( counts.shape[0], counts.shape[1], 1 ) )
    
    alpha = np.swapaxes( alpha, 0, 1 )
    alpha = alpha.reshape( ( 1, alpha.shape[0], alpha.shape[1] ) )
    
    alpha_sum = alpha.sum( axis=1 )
    counts_sum = counts.sum( axis=1 )
    
    p = gammaln( alpha_sum ) - gammaln( counts_sum + alpha_sum ) + \
        np.sum( gammaln( counts + alpha ), axis=1 ) - np.sum( gammaln( alpha ), axis=1 )
    
    return p
