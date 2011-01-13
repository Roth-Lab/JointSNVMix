'''
Created on 2010-11-14

@author: Andrew
'''
import numpy as np
from scipy.special import psi, polygamma
from scipy.optimize import fmin_ncg, fmin_bfgs, fmin_l_bfgs_b
from .log_pdf import log_beta_binomial_likelihood

def get_mle_p( vars ):
    x = vars[0]
    a = vars[1]
    b = vars[2]
    resp = vars[3]
    component = vars[4]
    
    return get_mle( x, a, b, resp, component )

def get_mle_2_p( vars ):
    x = vars[0]
    a = vars[1]
    b = vars[2]
    resp = vars[3]
    component = vars[4]
    
    return get_mle_2( x, a, b, resp, component )

def get_mom( a, b, resp, i ):
    d = a + b
    
    p_1 = ( resp * a ) / ( resp * d ) 
    p_1 = p_1[resp * d != 0]
    p_1_bar = np.mean( p_1 )
    p_1_sq_bar = np.mean( p_1 ** 2 )
    
    p_2 = ( resp * b ) / ( resp * d ) 
    p_2 = p_2[resp * d != 0]
    p_2_bar = np.mean( p_2 )
    
    s = ( p_1_bar - p_1_sq_bar ) / ( p_1_sq_bar - p_1_bar ** 2 )
    
    alpha = p_1_bar * s
    beta = p_2_bar * s
    
    if alpha <= 0:
        alpha = 1e-3
    if  beta <= 0:
        beta = 1e-3
    
    if i == 0:
        if alpha <= 1.:
            alpha = 1.01
        if beta >= 1.:
            beta = 0.99     
    if i == 1:
        if alpha <= 1.:
            alpha = 1.01
        if beta <= 1.:
            beta = 1.01
    if i == 2:
        if alpha >= 1.:
            alpha = 0.99
        if beta <= 1.:
            beta = 1.01

    x = np.array( [alpha, beta] )
    
    return x

def get_mle_fixed_point( x, a, b, resp, i ):
    alpha = x[0]
    beta = x[1]
    
    f = lambda y:-1 * get_f( y, 1, a, b, resp, i )

    f_new = f( x )
    diff = float( 'inf' )
    
    while diff > 1e-6:
        alpha_new = alpha * np.sum( resp * digamma_difference( a, alpha ) ) / np.sum( resp * digamma_difference( a + b, alpha + beta ) )
        beta_new = beta * np.sum( resp * digamma_difference( b, beta ) ) / np.sum( resp * digamma_difference( a + b, alpha + beta ) )
        
        alpha = alpha_new
        beta = beta_new
        
        x = np.array( [alpha, beta] )
        
        f_old = f_new
        f_new = f( x )
        
        diff = abs( f_new - f_old )
        
        print x
        print diff
    
    return x
        
def get_mle( x, a, b, resp, i ):
    if np.all( resp == 0 ):
        print "Empty class."
        return x
    
    if i == 0:
        bounds = [( 1, None ), ( 0, 1 )]
    elif i == 1:
        bounds = [( 1, None ), ( 1, None )]
    elif i == 2:
        bounds = [( 0, 1 ), ( 1, None )]
    
    f = lambda y:-1 * get_f( y, a, b, resp )
    g = lambda y:-1 * get_gradient( y, a, b, resp )
        
    x, r, s = fmin_l_bfgs_b( f, x, fprime=g, bounds=bounds )
    
    return x

def get_mle_2( x, a, b, resp, i ):
    if np.all( resp == 0 ):
        print "Empty class."
        return x
    
    t = 1.
    tolerance = 1e-6
    
    while True:
        print t
        f = lambda y:-1 * get_constrained_f( x, t, a, b, resp, i )
        g = lambda y:-1 * get_constrained_gradiet( x, t, a, b, resp, i )
        h = lambda y:-1 * get_constrained_hessian( x, t, a, b, resp, i )
        
#        x = fmin_ncg( f, x, g, fhess=h )
        x = fmin_newton( f, x, g, h )
        
        if 2 / t < tolerance:
            return x
        
        t = 20 * t

    return x

def fmin_newton( f, x, g, h, maxiter=80, tolerance=1e-3 ):
    iters = 0      

    while True:
        gx = g( x )
        hx = h( x )
        
        direction = np.linalg.solve( hx, gx )
        
        x_new, progress = newton_line_search( f, x, g, direction )
        
        if not progress:
            return x
        else:
            f_old = f( x )
            f_new = f( x_new ) 
            
            x = x_new
            
        diff = abs( f_new - f_old )
        print x, diff
        
        if iters > maxiter or diff < tolerance:
            return x
        
        iters += 1

def newton_line_search( f, x, g, direction, maxiter=20 ):
    step_size = 1.
    
    x_old = x
    x_new = x_old - step_size * direction
    
    f_old = f( x_old )
    g_old = g( x_old )
    
    f_new = f( x_new )
    
    iters = 0
    
    while f_old + 0.3 * step_size * np.dot( direction, g_old ) < f_new and iters < maxiter:
        step_size = 1 / 2. * step_size
        
        x_new = x_old - step_size * direction
        
        f_new = f( x_new )
        
        iters += 1
    
    if np.all( x_old == x_new ):
        print "No progress made."
        return x_old, False
    else:
        return x_new, True

def get_constrained_f( x, t, a, b, resp, i ):
    alpha = x[0]
    beta = x[1]
    
    if i == 0:
        phi = -np.log( ( alpha - 1. ) ) - np.log( beta )
    elif i == 1:
        phi = -np.log( ( alpha - 1 ) ) - np.log( ( beta - 1 ) )
    elif i == 2:
        phi = -np.log( alpha ) - np.log( ( beta - 1 ) )
        
    f_val = get_f( x, a, b, resp )
    
    f_val = t * f_val + phi
    
    return f_val

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

def get_constrained_gradiet( x, t, a, b, resp, i ):
    alpha = x[0]
    beta = x[1]
    
    if i == 0:
        grad_phi = np.array( [
                              1. / ( alpha - 1. ),
                              1. / beta
                              ] )
    elif i == 1:
        grad_phi = np.array( [
                              1. / ( alpha - 1. ),
                              1. / ( beta - 1. )
                              ] )
    elif i == 2:
        grad_phi = np.array( [
                              1. / ( alpha ),
                              1. / ( beta - 1. )
                              ] )
        
    grad = get_gradient( x, a, b, resp )
    
    grad = t * grad + grad_phi
    
    return grad


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

def get_constrained_hessian( x, t, a, b, resp, i ):
    alpha = x[0]
    beta = x[1]
    
    if i == 0:
        hess_phi = np.array( [
                              [ 1. / ( alpha - 1. ) ** 2, 0 ],
                              [ 0, 1. / beta ** 2]
                              ] )
    elif i == 1:
        hess_phi = np.array( [
                              [ 1. / ( alpha - 1. ) ** 2, 0 ],
                              [ 0, 1. / ( beta - 1. ) ** 2]
                              ] )
    elif i == 2:
        hess_phi = np.array( [
                              [ 1. / alpha ** 2, 0 ],
                              [ 0, 1. / ( beta - 1. ) ** 2]
                              ] )
        
    hess = get_hessian( x, a, b, resp )
    
    hess = t * hess + hess_phi
    
    return hess

def get_hessian( x, a, b, resp ):
    alpha = x[0]
    beta = x[1]
    
    d = a + b
    
    common_term = trigamma_difference( d, alpha + beta )
    
    deriv_wrt_alpha_and_alpha = trigamma_difference( a, alpha ) - common_term
    deriv_wrt_alpha_and_alpha = resp * deriv_wrt_alpha_and_alpha
    deriv_wrt_alpha_and_alpha = deriv_wrt_alpha_and_alpha.sum( axis=0 )
    
    deriv_wrt_beta_and_beta = trigamma_difference( b, beta ) - common_term
    deriv_wrt_beta_and_beta = resp * deriv_wrt_beta_and_beta
    deriv_wrt_beta_and_beta = deriv_wrt_beta_and_beta.sum( axis=0 )
    
    deriv_wrt_alpha_and_beta = -common_term
    deriv_wrt_alpha_and_beta = resp * deriv_wrt_alpha_and_beta
    deriv_wrt_alpha_and_beta = deriv_wrt_alpha_and_beta.sum( axis=0 )
    
    hess = np.array( [
                     [deriv_wrt_alpha_and_alpha, deriv_wrt_alpha_and_beta],
                     [deriv_wrt_alpha_and_beta, deriv_wrt_beta_and_beta]
                     ] )
    
    return hess
    
def digamma_difference( counts, parameter ):
    return psi( counts + parameter ) - psi( parameter )

def trigamma_difference( counts, parameter ):
    return polygamma( 1, counts + parameter ) - polygamma( 1, parameter )
