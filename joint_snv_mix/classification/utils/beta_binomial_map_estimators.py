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
    component = vars[4]
    alpha_prior = vars[5]
    beta_prior = vars[6]  
    
    return get_ml_estimates( x, a, b, resp, component, alpha_prior, beta_prior )

def get_ml_estimates( x, a, b, resp, i, alpha_prior, beta_prior ):
    if np.all( resp == 0 ):
        print "Empty class."
        return x
    
    f = lambda y:-1 * ( get_f( y, a, b, resp ) + \
                        get_penalty( y, alpha_prior, beta_prior ) )
    
    g = lambda y:-1 * ( get_gradient( y, a, b, resp ) + \
                        get_penalty_gradient( y, alpha_prior, beta_prior ) )
    
    # Bounds to keep the alpha, beta search away from 0.
    bounds = [( 1e-6, None ), ( 1e-6, None )]
    
    x = fmin_l_bfgs_b( f, x, fprime=g, bounds=bounds )
    
    return x[0]

def get_penalty( x, alpha_prior, beta_prior ):
    alpha = x[0]
    beta = x[1]
    
    alpha_prior_shape = alpha_prior[0]
    alpha_prior_scale = alpha_prior[1]
    alpha_prior_min = alpha_prior[2]
    
    beta_prior_shape = beta_prior[0]
    beta_prior_scale = beta_prior[1]
    beta_prior_min = beta_prior[2]
    
    alpha = alpha - alpha_prior_min
    beta = beta - beta_prior_min
    
    if alpha <= 0 or beta <= 0:
        return float( '-inf' )
    
    alpha_penalty = ( alpha_prior_shape - 1 ) * np.log( alpha ) - alpha / alpha_prior_scale
    beta_penalty = ( beta_prior_shape - 1 ) * np.log( beta ) - beta / beta_prior_scale
    
    return alpha_penalty + beta_penalty

def get_penalty_gradient( x, alpha_prior, beta_prior ):
    alpha = x[0]
    beta = x[1]
    
    alpha_prior_shape = alpha_prior[0]
    alpha_prior_scale = alpha_prior[1]
    alpha_prior_min = alpha_prior[2]
    
    beta_prior_shape = beta_prior[0]
    beta_prior_scale = beta_prior[1]
    beta_prior_min = beta_prior[2]
    
    alpha = alpha - alpha_prior_min
    beta = beta - beta_prior_min

    alpha_penalty = ( alpha_prior_shape - 1 ) / alpha - 1 / alpha_prior_scale
    beta_penalty = ( beta_prior_shape - 1 ) / beta - 1 / beta_prior_scale
    
    grad_penalty = np.array( [alpha_penalty, beta_penalty] )
    
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
