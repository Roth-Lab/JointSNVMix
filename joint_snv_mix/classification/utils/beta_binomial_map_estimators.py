'''
Created on 2010-11-14

@author: Andrew
'''
import numpy as np
from scipy.special import psi
from scipy.optimize import fmin_l_bfgs_b
from .log_pdf import log_beta_binomial_likelihood

def get_mle_p( vars ):
    x = vars[0]
    a = vars[1]
    b = vars[2]
    resp = vars[3]
    alpha_prior = vars[4]
    beta_prior = vars[5]  
    
    return get_ml_estimates( x, a, b, resp, alpha_prior, beta_prior )

def get_ml_estimates( x, a, b, resp, precision_prior, location_prior ):
    if np.all( resp == 0 ):
        print "Empty class."
        return x
    
    f = lambda y:-1 * ( get_f( y, a, b, resp ) + \
                        get_penalty( y, precision_prior, location_prior ) )
    
    g = lambda y:-1 * ( get_gradient( y, a, b, resp ) + \
                        get_penalty_gradient( y, precision_prior, location_prior ) )
    
    # Bounds to keep the alpha, beta search away from 0.
    bounds = [( 1e-6, None ), ( 1e-6, None )]
    
    x = fmin_l_bfgs_b( f, x, fprime=g, bounds=bounds )
    
    return x[0]

def get_penalty( x, precision_prior, location_prior ):
    alpha = x[0]
    beta = x[1]
    
    precision_prior_shape = precision_prior[0]
    precision_prior_scale = precision_prior[1]
    precision_prior_min = precision_prior[2]
    
    location_prior_alpha = location_prior[0]
    location_prior_beta = location_prior[1]

    s = alpha + beta
    mu = alpha / s
    
    s = s - precision_prior_min
    
    if s <= 0:
        return float( '-inf' )
    
    scale_penalty = ( precision_prior_shape - 1 ) * np.log( s ) - s / precision_prior_scale
    location_penalty = ( location_prior_alpha - 1 ) * np.log( mu ) + ( location_prior_beta - 1 ) * np.log( 1 - mu )
    
    return scale_penalty + location_penalty

def get_penalty_gradient( x, precision_prior, location_prior ):
    alpha = x[0]
    beta = x[1]
    
    precision_prior_shape = precision_prior[0]
    precision_prior_scale = precision_prior[1]
    precision_prior_min = precision_prior[2]
    
    location_prior_alpha = location_prior[0]
    location_prior_beta = location_prior[1]

    s = alpha + beta
    mu = alpha / s
    
    s = s - precision_prior_min

    grad_scale_penalty = ( precision_prior_shape - 1 ) / s - 1 / precision_prior_scale
    grad_location_penalty = ( location_prior_alpha - 1 ) / mu - ( location_prior_beta - 1 ) / ( 1 - mu )
    
    grad_penalty = np.array( [grad_scale_penalty, grad_location_penalty] )
    
    return grad_penalty

def get_f( x, a, b, resp ):
    alpha = x[0]
    beta = x[1]
    
    nrows = a.shape[0]
    
    d = a + b
    
    f_val = log_beta_binomial_likelihood( a, d, alpha, beta )
    f_val = f_val.reshape( ( nrows, ) )
    f_val = resp * f_val
    f_val = f_val.sum( axis=0 ) 
    
    return f_val

def get_gradient( x, a, b, resp ):
    alpha = x[0]
    beta = x[1]
    
    d = a + b
    
    common_term = digamma_difference( d, alpha + beta )
    
    deriv_wrt_alpha = digamma_difference( a, alpha ) - common_term
    deriv_wrt_alpha = resp * deriv_wrt_alpha
    deriv_wrt_alpha = deriv_wrt_alpha.sum( axis=0 )
    
    deriv_wrt_beta = digamma_difference( b, beta ) - common_term
    deriv_wrt_beta = resp * deriv_wrt_beta
    deriv_wrt_beta = deriv_wrt_beta.sum( axis=0 )
    
    grad = np.array( [deriv_wrt_alpha, deriv_wrt_beta] ) 
    
    return grad
    
def digamma_difference( counts, parameter ):
    return psi( counts + parameter ) - psi( parameter )
