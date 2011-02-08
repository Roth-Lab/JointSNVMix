'''
Created on 2010-11-14

@author: Andrew
'''
import numpy as np
from scipy.special import psi
from scipy.optimize import fmin_l_bfgs_b
from .log_pdf import log_beta_binomial_likelihood

def get_junk_map( data, params, priors, resp ):        
    prec_shape = priors['precision']['shape']
    prec_scale = priors['precision']['scale']
    prec_min = priors['precision']['min']
    
    loc_alpha = priors['location']['alpha']
    loc_beta = priors['location']['beta']
    
        
    f = lambda y:-1 * ( get_junk_f( y, data, resp ) + \
                        get_penalty( y, loc_alpha, loc_beta, prec_shape, prec_scale, prec_min ) )
    
    g = lambda y:-1 * ( get_junk_gradient( y, data, resp ) + \
                        get_penalty_gradient( y, loc_alpha, loc_beta, prec_shape, prec_scale, prec_min ) )
    
    # Bounds to keep the alpha, beta search away from 0.
    bounds = [( 1e-6, None ), ( 1e-6, None )]
    
    x = ( params['alpha'], params['alpha'] )
    
    x = fmin_l_bfgs_b( f, x, fprime=g, bounds=bounds )
    
    params = {}
    params['alpha'] = x[0][0]
    params['beta'] = x[0][1]
    
    return params
    
def get_junk_f( x, data, resp ):
    normal_f = get_f( data.a['normal'], data.b['normal'], resp )
    tumour_f = get_f( data.a['tumour'], data.b['tumour'], resp )
    
    return normal_f + tumour_f

def get_junk_gradient( x, data, resp ):
    normal_f = get_gradient( data.a['normal'], data.b['normal'], resp )
    tumour_f = get_gradient( data.a['tumour'], data.b['tumour'], resp )
    
    return normal_f + tumour_f



def get_mle_p( vars ):
    x = vars[0]
    a = vars[1]
    b = vars[2]
    resp = vars[3]
    location_prior = vars[4]
    precision_prior = vars[5]
    component = vars[6]
    
    return get_ml_estimates( x, a, b, resp, location_prior, precision_prior, component )

def get_ml_estimates( x, a, b, resp, location_prior, precision_prior, component ):
    if np.all( resp == 0 ):
        print "Empty class."
        return x
    
    prec_shape = precision_prior['shape'][component]
    prec_scale = precision_prior['scale'][component]
    prec_min = precision_prior['min'][component]
    
    loc_alpha = location_prior['alpha'][component]
    loc_beta = location_prior['beta'][component]
    
    f = lambda y:-1 * ( get_f( y, a, b, resp ) + \
                        get_penalty( y, loc_alpha, loc_beta, prec_shape, prec_scale, prec_min ) )
    
    g = lambda y:-1 * ( get_gradient( y, a, b, resp ) + \
                        get_penalty_gradient( y, loc_alpha, loc_beta, prec_shape, prec_scale, prec_min ) )
    
    # Bounds to keep the alpha, beta search away from 0.
    bounds = [( 1e-6, None ), ( 1e-6, None )]
    
    x = fmin_l_bfgs_b( f, x, fprime=g, bounds=bounds )
    
    return x[0]

def get_penalty( x, loc_alpha, loc_beta, prec_shape, prec_scale, prec_min ):
    alpha = x[0]
    beta = x[1]
    
    mu = alpha / ( alpha + beta )
    s = alpha + beta
    
    s = s - prec_min
    
    if s <= 0:
        return float( '-inf' )
    
    scale_penalty = ( prec_shape - 1 ) * np.log( s ) - s / prec_scale
    location_penalty = ( loc_alpha - 1 ) * np.log( mu ) + ( loc_beta - 1 ) * np.log( 1 - mu )
    
    return scale_penalty + location_penalty

def get_penalty_gradient( x, loc_alpha, loc_beta, prec_shape, prec_scale, prec_min ):
    alpha = x[0]
    beta = x[1]
    
    mu = alpha / ( alpha + beta )
    s = alpha + beta
    
    s = s - prec_min

    d_pen_s = ( prec_shape - 1 ) / s - 1 / prec_scale
    d_pen_mu = ( loc_alpha - 1 ) / mu - ( loc_beta - 1 ) / ( 1 - mu )
    
    d_mu_alpha = beta / ( alpha + beta ) ** 2
    d_mu_beta = -alpha / ( alpha + beta ) ** 2
    
    d_pen_alpha = d_pen_mu * d_mu_alpha + d_pen_s
    d_pen_beta = d_pen_mu * d_mu_beta + d_pen_s
    
    grad_penalty = np.array( [d_pen_alpha, d_pen_beta] )
    
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
