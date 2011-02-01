'''
Created on 2011-01-31

@author: Andrew Roth
'''
import numpy as np

from scipy.special import betaln, psi
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize.tnc import fmin_tnc
from joint_snv_mix import constants
import ConfigParser
from scipy.cluster.vq import kmeans2
from scipy.optimize.optimize import check_grad


class EMTrainer( object ):
    def __init__( self, tolerance=1e-6, max_iters=1000 ):
        self.tolerance = tolerance
        self.max_iters = max_iters
        
    def train( self, density, data, resp ):
        iters = 0
        converged = False
        
        old_posterior_value = float( '-inf' )
  
        while not converged:            
            density.update( data, resp )
            
            resp = density.get_responsibilities( data )

            posterior_value = density.log_likelihood( data )
            
            if iters > 0:
                posterior_change = ( posterior_value - old_posterior_value ) / abs( old_posterior_value )
            else:
                posterior_change = float( 'inf' )
            
            self._print_diagnostic_message( iters, posterior_value, old_posterior_value, posterior_change )
            old_posterior_value = posterior_value
            
            if posterior_change < 0:
                print "Posterior decreased. This could be a bug or overly stringent convergence criterion."
                converged = True
            
            elif posterior_change < self.tolerance:
                converged = True
                            
            if iters >= self.max_iters:
                print "Maximum numbers of EM iterations exceeded. Exiting training."                
                converged = True
            
            iters += 1
                     
        return density
    
    def _print_diagnostic_message( self, iters, posterior_value, old_posterior_value, posterior_change ):
        print "#" * 100
        print "# Diagnostics."
        print "#" * 100
        print "Number of iterations : ", iters
        print "New posterior : ", posterior_value
        print "Old posterior : ", old_posterior_value 
        print "Posterior change : ", posterior_change
    
        print "Parameters :"
        
#        for param_name, param_value in self.posterior.parameters.items():
#            print param_name, param_value         

def initialise_joint_responsibilities( data ):

    '''
    Intialise responsibilities via k-means clustering.
    '''
    a_1 = np.asarray( data.a[0], dtype=np.float64 )
    b_1 = np.asarray( data.b[0], dtype=np.float64 )
    p_1 = a_1 / ( a_1 + b_1 )
          
    a_2 = np.asarray( data.a[1], dtype=np.float64 )
    b_2 = np.asarray( data.b[1], dtype=np.float64 )
    p_2 = a_2 / ( a_2 + b_2 )

    shape = ( data.nrows, 9 )
    
    responsibilities = np.zeros( shape )
    
    init_centers = np.array( ( 1., 0.5, 0. ) )
    
    cluster_centers_1, labels_1 = kmeans2( p_1, init_centers, minit='matrix' )
    cluster_centers_2, labels_2 = kmeans2( p_2, init_centers, minit='matrix' )

    labels = 3 * labels_1 + labels_2

    for id in range( 9 ):
        index = labels == id
        
        responsibilities[index, id] = 1.
    
    return responsibilities

#=======================================================================================================================
# Beta-Binomial Code
#=======================================================================================================================
class SixParameterBetaBinomialMixture( object ):
    def __init__( self, priors_file=None, params_file=None ):
        if priors_file is not None:
            self._init_from_priors( priors_file )
            
            self.can_train = True
        else:
            self._init_from_parameters( params_file )
            
            self.can_train = False

    def train( self, data ):
        if not self.can_train:
            raise Exception( "Can't train without specifiying a priors file." )
        
        init_resp = initialise_joint_responsibilities( data )
        
        trainer = EMTrainer()
        
        self.mixture_density = trainer.train( self.mixture_density, data, init_resp )
        
    def classify( self, data ):
        return self.mixture_density.get_responsibilities( data )

    def _init_from_priors( self, priors_file_name ):
        priors = self._parse_priors_file( priors_file_name ) 
        
        density_priors = {} 
        
        density_priors['normal'] = BetaBinomialPrior( priors['densities']['normal'] )
        
        density_priors['tumour'] = BetaBinomialPrior( priors['densities']['tumour'] )

        density_params = { 'normal' : {}, 'tumour' : {} }
                
        density_params['normal']['alpha'] = np.array( [99, 50, 1], np.float )
        density_params['tumour']['alpha'] = np.array( [99, 50, 1], np.float )
        
        density_params['normal']['beta'] = np.array( [1, 50, 99], np.float )
        density_params['tumour']['beta'] = np.array( [1, 50, 99], np.float )
        
        component_densities = JointBetaBinomialDensities( density_priors, density_params )
        
        mw_priors = MixWeightsPriors( priors['mix-weights'] )
        
        mix_weights = MixWeights( mw_priors ) 
        
        self.mixture_density = JointBetaBinomialMixtureDensity( component_densities, mix_weights )
        
    def _init_from_parameters( self, parameters_file_name ):
        parameters = self._parse_parameters_file( parameters_file_name )
        
        priors = {'normal' : None, 'tumour' : None}
        mw_priors = None
        
        density_parameters = parameters['densities']
        mix_weight_value = parameters['mix-weights']
        
        component_densities = JointBetaBinomialDensities( priors, density_parameters )
        
        mix_weights = MixWeights( mw_priors, init_value=mix_weight_value ) 
        
        self.mixture_density = JointBetaBinomialMixtureDensity( component_densities, mix_weights )
        
    def _parse_priors_file( self, priors_file_name ):
        parser = ConfigParser.SafeConfigParser()
        parser.read( priors_file_name )
        
        priors = {}
        
        priors['mix-weights'] = np.zeros( ( 9, ) )
    
        for i, genotype_tuple in enumerate( constants.joint_genotypes ):
            genotype = "_".join( genotype_tuple )
        
            priors['mix-weights'][i] = parser.getfloat( 'mix-weights', genotype )
        
        priors['densities'] = {}        
        for genome in constants.genomes:
            priors['densities'][genome] = {}
            
            priors['densities'][genome]['location'] = {}
            priors['densities'][genome]['location']['alpha'] = np.zeros( ( 3, ) )
            priors['densities'][genome]['location']['beta'] = np.zeros( ( 3, ) )
        
            priors['densities'][genome]['precision'] = {}
            priors['densities'][genome]['precision']['shape'] = np.zeros( ( 3, ) )
            priors['densities'][genome]['precision']['scale'] = np.zeros( ( 3, ) )
            priors['densities'][genome]['precision']['min'] = np.zeros( ( 3, ) )
            
            for i, genotype in enumerate( constants.genotypes ):                                
                genome_genotype = "_".join( ( genome, genotype ) )
            
                priors['densities'][genome]['location']['alpha'][i] = parser.getfloat( 'location_alpha', genome_genotype )
                priors['densities'][genome]['location']['beta'][i] = parser.getfloat( 'location_beta', genome_genotype )
                
                priors['densities'][genome]['precision']['shape'][i] = parser.getfloat( 'precision_shape', genome_genotype )
                priors['densities'][genome]['precision']['scale'][i] = parser.getfloat( 'precision_scale', genome_genotype )
                priors['densities'][genome]['precision']['min'][i] = parser.getfloat( 'precision_min', genome_genotype )
            
        
        return priors
    
    def _parse_parameters_file( self, parameters_file_name ):
        parser = ConfigParser.SafeConfigParser()
        parser.read( parameters_file_name )
        
        parameters = {}
                
        parameters['mix-weights'] = np.zeros( ( 9, ) )
    
        for i, genotype_tuple in enumerate( constants.joint_genotypes ):
            genotype = "_".join( genotype_tuple )
        
            parameters['mix-weights'][i] = parser.getfloat( 'mix-weights', genotype )

        parameters['densities'] = {}
        for genome in constants.genomes:
            parameters['densities'][genome] = {}
            parameters['densities'][genome]['alpha'] = np.zeros( ( 3, ) )
            parameters['densities'][genome]['beta'] = np.zeros( ( 3, ) )
                        
            for i, genotype in enumerate( constants.genotypes ):
                genome_genotype = "_".join( ( genome, genotype ) )
            
                parameters['densities'][genome]['alpha'][i] = parser.getfloat( 'alpha', genome_genotype )
                parameters['densities'][genome]['beta'][i] = parser.getfloat( 'beta', genome_genotype )
            
        return parameters

class JointBetaBinomialMixtureDensity( object ):
    def __init__( self, component_densities, mix_weights, ):
        self.mix_weights = mix_weights
        
        self.densities = component_densities
    
    def get_responsibilities( self, data ):
        log_resp = self.densities.log_likelihood( data )
        
        log_resp = log_resp + np.log( self.mix_weights.value )
        
        log_resp = self._log_space_normalise_rows( log_resp )
        
        return np.exp( log_resp )
        
    def update( self, data, responsibilities ):
        self.mix_weights.update( responsibilities )
        
        self.densities.update( data, responsibilities )
        
    def log_likelihood( self, data ):
        ll = 0
        
        ll += self.mix_weights.priors_log_likelihood()
        
        ll += np.sum( self.densities.priors_log_likelihood() )
        
        ll += np.sum( self.densities.log_likelihood( data ) + np.log( self.mix_weights.value ) )
        
        return ll
    
    def _log_space_normalise_rows( self, log_X ):
        nrows = log_X.shape[0]
        shape = ( nrows, 1 )
        
        log_norm_const = np.logaddexp.reduce( log_X, axis=1 )
        log_norm_const = log_norm_const.reshape( shape )
    
        log_X = log_X - log_norm_const
                
        return log_X

class JointBetaBinomialDensities( object ):
    def __init__( self, priors, init_params ):       
        self.n_normal_class = init_params['normal']['alpha'].size
        self.n_tumour_class = init_params['tumour']['alpha'].size
        
        self.normal_compents = BetaBinomialMixtureVariable( init_params['normal'], priors['normal'] )
        
        self.tumour_compents = BetaBinomialMixtureVariable( init_params['tumour'], priors['tumour'] )
    
    def log_likelihood( self, data ):
        normal_ll = self.normal_compents.log_likelihood( data.a[0], data.b[0] )
        tumour_ll = self.tumour_compents.log_likelihood( data.a[1], data.b[1] )
        
        shape = ( data.nrows, 1 )
        
        ll = [ 
              normal_ll[:, i].reshape( shape ) + tumour_ll for i in range( self.n_normal_class )
             ]
        
        ll = np.hstack( ll ) 
        
        return ll
    
    def priors_log_likelihood( self ):
        p_ll = 0
        
        p_ll += self.normal_compents.priors_log_likelihood()
        p_ll += self.tumour_compents.priors_log_likelihood()
        
        return p_ll
    
    def update( self, data, responsibilities ):
        normal_resp, tumour_resp = self._get_marginals( responsibilities )
        
        self.normal_compents.update( data.a[0], data.b[0], normal_resp )
        self.tumour_compents.update( data.a[1], data.b[1], tumour_resp )
        
    def _get_marginals( self, responsibilities ):
        nrows = responsibilities.shape[0]
    
        shape = ( nrows, self.n_normal_class, self.n_tumour_class )
        
        responsibilities = responsibilities.reshape( shape )
        
        normal_marginals = responsibilities.sum( axis=2 )
        tumour_marginals = responsibilities.sum( axis=1 )
        
        return normal_marginals, tumour_marginals

class BetaBinomialMixtureVariable( object ):
    def __init__( self, params, priors=None ):
        alpha = params['alpha']
        beta = params['beta']
        
        if np.isscalar( alpha ):        
            alpha = np.array( alpha )
            beta = np.array( beta )
        
        self.priors = priors
        
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
    
    def priors_log_likelihood( self ):
        if self.priors is None:
            return 0
        
        return self.priors.log_likelihood( self.alpha, self.beta )
    
    def update( self, a, b, resp, penalty_func=None ):        
        vars = []
        
        alpha = self.alpha.flatten()
        beta = self.beta.flatten()
        
        for i in range( self.nclass ):
            x0 = ( alpha[i], beta[i] )
            
            r = resp[:, i]
        
            f_unpen = lambda x:-1 * self._complete_data_log_likelihood( a, b, r, alpha=x[0], beta=x[1] )
            g_unpen = lambda x:-1 * self._complete_data_log_likelihood_gradient( a, b, r, alpha=x[0], beta=x[1] )
            
            if self.use_penalty:
                f = lambda x: f_unpen( x ) - self.penalty_func( x[0], x[1], component=i )
                g = lambda x: g_unpen( x ) - self.penalty_grad( x[0], x[1], component=i )
            else:
                f = f_unpen
                g = g_unpen
            
            print check_grad( f, g, x0 )
            
            vars.append( 
                        ( x0, f, g )
                        )

        self._run_update( vars )
    
    def _run_update( self, vars ):
        alpha = np.zeros( ( 3, ) )
        beta = np.zeros( ( 3, ) )
        
        bounds = [( 1e-6, None ), ( 1e-6, None )]
        
        for i, var in enumerate( vars ):
            x = fmin_l_bfgs_b( var[1], var[0], fprime=var[2], bounds=bounds )
        
            alpha[i] = x[0][0]
            beta[i] = x[0][1]
        
            print x
        
        self.alpha = self._reshape_parameters( alpha )
        self.beta = self._reshape_parameters( beta )
    
    def _complete_data_log_likelihood( self, a, b, resp, alpha=None, beta=None ):     
        if alpha <= 0 or beta <= 0:
            return float( "-inf" )
           
        f_val = self.log_likelihood( a, b, alpha=alpha, beta=alpha )
        f_val = resp.flatten() * f_val.flatten()
        f_val = f_val.sum()

        return f_val
    
    def _complete_data_log_likelihood_gradient( self, a, b, resp, alpha=None, beta=None ):        
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
    def __init__( self, hyper_parameters ):       
        self.location = hyper_parameters['location']
        self.precision = hyper_parameters['precision']
    
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
    
#=======================================================================================================================
# Mix-weights
#=======================================================================================================================
class MixWeights( object ):
    def __init__( self, priors, init_value=None ):
        self.priors = priors
        
        self.value = init_value
        
    def update( self, responsibilities ):
        N_g = responsibilities.sum( axis=0 )

        mix_weights = N_g + self.priors.value - 1

        mix_weights = np.exp( np.log( mix_weights ) - np.log( mix_weights.sum() ) )

        self.value = mix_weights
    
    def priors_log_likelihood( self ):
        return self.priors.log_likelihood( self.value )
        
class MixWeightsPriors( object ):
    def __init__( self, value ):
        self.value = value
        
    def log_likelihood( self, mix_weight ):                
        ll = ( self.value - 1 ) * np.log( mix_weight )
        
        ll = np.sum( ll )
        
        return ll
