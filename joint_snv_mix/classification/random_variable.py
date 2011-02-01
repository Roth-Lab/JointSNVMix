'''
Created on 2011-01-31

@author: Andrew Roth
'''
import numpy as np

from scipy.special import betaln, psi
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize.tnc import fmin_tnc

class JointBetaBinomialMixtureDensity( object ):
    def __init__( self, mix_weights, priors, parameters ):
#        self.mix_weight_priors = pass
        
#        self.mix_weights = MixWeights( mix_weights, priors['mix_weights'] )
        
        self.densities = JointBetaBinomialDensities( priors['densities']['normal'], priors['densities']['tumour'] )
        
    def update( self, data, responsibilities ):
        self.mix_weights.update( data, responsibilities )
        
        self.densities.update( data, responsibilities )
        
    def log_likelihood( self, data ):
        ll = 0
        
        ll += self.mix_weights.priors.log_likelihood()
        
        ll += np.sum( self.densities.normal_priors.log_likelihood() )
        ll += np.sum( self.densities.tumour_priors.log_likelihood() )
        
        ll += np.sum( self.densities.log_likelihood() + np.log( self.mix_weights.value ) )
        
        return ll

class JointBetaBinomialDensities( object ):
    def __init__( self, normal_priors, tumour_priors ):
        self.n_normal_class = 3
        self.n_tumour_class = 3
        
        self.normal_priors = normal_priors
        self.tumour_priors = tumour_priors
        
        alpha = np.array( [99, 5, 1], np.float )
        beta = np.array( [1, 5, 99], np.float )
        
        self.normal_compents = BetaBinomialMixtureVariable( alpha, beta, normal_priors )
        self.tumour_compents = BetaBinomialMixtureVariable( alpha, beta, tumour_priors )
    
    def log_likelihood( self, data ):
        normal_ll = self.normal_compents.log_likelihood( data.a[0], data.b[0] )
        tumour_ll = self.tumour_compents.log_likelihood( data.a[1], data.b[1] )
        
        shape = ( data.nrows, 1 )
        
        ll = [ 
              normal_ll[:, i].reshape( shape ) + tumour_ll for i in range( self.n_normal_class )
             ]
        
        ll = np.hstack( ll ) 
        
        return ll
    
    def update( self, data, responsibilities ):
        normal_resp, tumour_resp = self._get_marginals( responsibilities )
        
        self.normal_compents.update( data.a[0], data.b[0], normal_resp )
        self.normal_compents.update( data.a[1], data.b[1], tumour_resp )
        
    def _get_marginals( self, responsibilities ):
        nrows = responsibilities.shape[0]
    
        shape = ( nrows, self.n_normal_class, self.n_tumour_class )
        
        responsibilities = responsibilities.reshape( shape )
        
        normal_marginals = responsibilities.sum( axis=2 )
        tumour_marginals = responsibilities.sum( axis=1 )
        
        return normal_marginals, tumour_marginals


class BetaBinomialMixtureVariable( object ):
    def __init__( self, alpha, beta, priors=None ):
        if np.isscalar( alpha ):        
            alpha = np.array( alpha )
            beta = np.array( beta )
        
        if priors is None:
            self.use_penalty = False
        else:
            self.use_penalty = True
            
            self.penalty_func = priors.log_likelihood
            self.penalty_grad = priors.log_likelihood_grad
        
        self.nclass = alpha.size

        self.alpha = self._reshape_parameters( alpha )
        self.beta = self._reshape_parameters( beta )
    
    def log_likelihood( self, a, b, alpha=None, beta=None ):
        a = self._reshape_counts( a )
        b = self._reshape_counts( b )
        
        if alpha is None:
            alpha = self.alpha
        if beta is None:
            beta = self.beta
        
        ll = betaln( a + alpha, b + beta ) - betaln( alpha, beta )
        
        return ll
    
    def update( self, a, b, resp, penalty_func=None ):        
        vars = []
        
        alpha = self.alpha.flatten()
        beta = self.beta.flatten()
        
        for i in range( self.nclass ):
            x0 = ( alpha[i], beta[i] )
            
            r = resp[:, i]
        
            f_unpen = lambda x:-1 * self._mixture_log_likelihood( a, b, r, alpha=x[0], beta=x[1] )
            g_unpen = lambda x:-1 * self._mixture_log_likelihood_gradient( a, b, r, alpha=x[0], beta=x[1] )
            
            if self.use_penalty:
                f = lambda x: f_unpen( x ) - self.penalty_func( x[0], x[1], component=i )
                g = lambda x: g_unpen( x ) - self.penalty_grad( x[0], x[1], component=i )
            else:
                f = f_unpen
                g = g_unpen
            
            vars.append( 
                        ( x0, f, g )
                        )

        self._run_update( vars )
    
    def _run_update( self, vars ):
        alpha = np.zeros( ( 3, ) )
        beta = np.zeros( ( 3, ) )
        
        bounds = [( 1e-6, None ), ( 1e-6, None )]
        
        for i, var in enumerate( vars ):
            x = fmin_tnc( var[1], var[0], fprime=var[2], bounds=bounds )
        
            alpha[i] = x[0][0]
            beta[i] = x[0][1]
        
            print x
        
        self.alpha = self._reshape_parameters( alpha )
        self.beta = self._reshape_parameters( beta )
    
    def _mixture_log_likelihood( self, a, b, resp, alpha=None, beta=None ):     
        if alpha <= 0 or beta <= 0:
            return float( "-inf" )
           
        f_val = self.log_likelihood( a, b, alpha=alpha, beta=alpha )
        f_val = resp.flatten() * f_val.flatten()
        f_val = f_val.sum()

        return f_val
    
    def _mixture_log_likelihood_gradient( self, a, b, resp, alpha=None, beta=None ):        
        d = a + b
        
        pochammer = self._digamma_difference
        
        common_term = pochammer( d, alpha + beta )
        
        deriv_wrt_alpha = pochammer( a, alpha ) - common_term
        deriv_wrt_alpha = resp * deriv_wrt_alpha
        deriv_wrt_alpha = deriv_wrt_alpha.sum( axis=0 )
        
        deriv_wrt_beta = pochammer( b, beta ) - common_term
        deriv_wrt_beta = resp * deriv_wrt_beta
        deriv_wrt_beta = deriv_wrt_beta.sum( axis=0 )
        
        grad = np.array( [deriv_wrt_alpha, deriv_wrt_beta] ) 
        
        return grad

    def _digamma_difference( self, counts, parameter ):
        return psi( counts + parameter ) - psi( parameter )
    
    def _reshape_counts( self, counts ):
        n = counts.size
        counts_shape = ( n, 1 )
        
        counts = counts.reshape( counts_shape )
        
        return counts
    
    def _reshape_parameters( self, parameters ):
        params_shape = ( 1, self.nclass )
        
        parameters = parameters.reshape( params_shape )
        
        return parameters
    
class BetaBinomialPrior( object ):
    def __init__( self, location, precision ):
        self.location = location
        self.precision = precision
    
    def log_likelihood( self, alpha, beta, component=None ):
        s = alpha + beta
        mu = alpha / s
        
        if component is None:
            precision_min = self.precision['min']
            precision_shape = self.precision['shape']
            precision_scale = self.precision['scale']        
            location_alpha = self.location['alpha']
            location_beta = self.location['beta']
        else:
            precision_min = self.precision['min'][component]
            precision_shape = self.precision['shape'][component]
            precision_scale = self.precision['scale'][component]
            location_alpha = self.location['alpha'][component]
            location_beta = self.location['beta'][component]
        
        s = s - precision_min
        
        if np.isscalar( s ) <= 0:
            return float( '-inf' )
        
        scale_penalty = ( precision_shape - 1 ) * np.log( s ) - s / precision_scale
        location_penalty = ( location_alpha - 1 ) * np.log( mu ) + ( location_beta - 1 ) * np.log( 1 - mu )
        
        return scale_penalty + location_penalty
    
    def log_likelihood_grad( self, alpha, beta, component=None ):
        s = alpha + beta
        mu = alpha / s
        
        if component is None:
            precision_min = self.precision['min']
            precision_shape = self.precision['shape']
            precision_scale = self.precision['scale']        
            location_alpha = self.location['alpha']
            location_beta = self.location['beta']
        else:
            precision_min = self.precision['min'][component]
            precision_shape = self.precision['shape'][component]
            precision_scale = self.precision['scale'][component]
            location_alpha = self.location['alpha'][component]
            location_beta = self.location['beta'][component]

        s = s - precision_min
    
        grad_scale_penalty = ( precision_shape - 1 ) / s - 1 / precision_scale
        grad_location_penalty = ( location_alpha - 1 ) / mu - ( location_beta - 1 ) / ( 1 - mu )
        
        grad_penalty = np.array( [grad_scale_penalty, grad_location_penalty] )
        
        return grad_penalty
