'''
Created on 2011-01-18

@author: Andrew Roth
'''
import numpy as np
import multiprocessing
from joint_snv_mix.classification.utils.beta_binomial_map_estimators import get_mle_p

def get_marginals( responsibilities, nclass ):
    nrows = responsibilities.shape[0]

    shape = ( nrows, nclass, nclass )
    
    responsibilities = responsibilities.reshape( shape )
    
    normal_marginals = responsibilities.sum( axis=2 )
    tumour_marginals = responsibilities.sum( axis=1 )
    
    marginals = []

    marginals.append( normal_marginals )
    marginals.append( tumour_marginals )

    return marginals

#=======================================================================================================================
# Abstract Models
#=======================================================================================================================
class EMPosterior( object ):
    def __init__( self, data, priors, responsibilities ):
        self.data = data
        self.priors = priors
        self.responsibilities = responsibilities
        
        self._init_parameters()
    
    def _init_parameters( self ):
        raise NotImplemented

    def update( self, responsibilities ):
        self.responsibilities = responsibilities     

        self._update_mix_weights()
        
        self._update_density_parameters()

    def _update_mix_weights( self ):
        N_g = self.responsibilities.sum( axis=0 )

        mix_weights = N_g + self.priors['kappa'] - 1

        mix_weights = np.exp( np.log( mix_weights ) - np.log( mix_weights.sum() ) )

        self.parameters['pi'] = mix_weights
        
    def _update_density_parameters( self ):
        raise NotImplemented
    
#=======================================================================================================================
# Independent Models
#=======================================================================================================================
class IndependentBetaBinomialPosterior( EMPosterior ):
    def __init__( self, data, priors, responsibilities ):
        EMPosterior.__init__( self, data, priors, responsibilities )
        
        self.pool = multiprocessing.Pool( maxtasksperchild=1 )
    
    def _init_parameters( self ):
        '''
        Initialise parameters. This is only necessary to initialise gradient descent.
        '''
        self.parameters = {}
        
        self._update_mix_weights()
        
        self.parameters['alpha'] = np.array( [99, 5, 1], np.float )

        self.parameters['beta'] = np.array( [1, 5, 99] , np.float )
        
        print "Initial parameter values : ", self.parameters
    
    def _update_density_parameters( self ):        
        a = self.data.a
        b = self.data.b
        
        vars = []
              
        for component in range( 3 ):
            x = np.zeros( ( 2, ) )
            
            x[0] = self.parameters['alpha'][component]
            x[1] = self.parameters['beta'][component]
            
            resp = self.responsibilities[:, component]
            
            precision_prior = self.priors['precision'][component]
            location_prior = self.priors['location'][component]
            
            vars.append( [x, a, b, resp, component, precision_prior, location_prior] )
                
        results = self.pool.map( get_mle_p, vars )
        
        for component in range( 3 ):
            self.parameters['alpha'][component] = results[component][0]
            self.parameters['beta'][component] = results[component][1]

class IndependentBinomialPosterior( EMPosterior ):
    def _init_parameters( self ):
        self.parameters = {}
        
        self._update_mix_weights()

        alpha = self.priors['alpha']
        beta = self.priors['beta']
        
        self.parameters['mu'] = alpha / ( alpha + beta )
        
        print self.parameters

    def _update_density_parameters( self ):
        a = self.data.a
        b = self.data.b
        d = a + b
        
        shape = ( self.data.nrows, 1 )
        a = a.reshape( shape )
        d = d.reshape( shape )
        
        alpha = self.priors['alpha']
        beta = self.priors['beta']
        
        resp = self.responsibilities
        
        a_bar = np.sum( resp * a, axis=0 )
        
        d_bar = np.sum( resp * d, axis=0 )
        
        numerator = a_bar + alpha - 1
        denominator = d_bar + alpha + beta - 2
        
        self.parameters['mu'] = numerator / denominator
        
#=======================================================================================================================
# Joint Models
#=======================================================================================================================
class JointBetaBinomialPosterior( EMPosterior ):
    def __init__( self, data, priors, responsibilities, nclass=3 ):
        self.nclass = nclass
        
        ncpus = nclass * 2
        
        self.pool = multiprocessing.Pool( processes=ncpus, maxtasksperchild=1 )
        
        EMPosterior.__init__( self, data, priors, responsibilities )
    
    def _init_parameters( self ):
        '''
        Initialise parameters. This is only necessary to initialise gradient descent. 
        '''
        self.parameters = {}
        
        self._update_mix_weights()

        location_alpha = self.priors['location'][:, :, 0]
        location_beta = self.priors['location'][:, :, 1]
        
        precision_shape = self.priors['precision'][:, :, 0]
        precision_scale = self.priors['precision'][:, :, 1] 
        
        s = precision_shape * precision_scale
        mu = location_alpha / ( location_alpha + location_beta ) 
        
        self.parameters['alpha'] = s * mu
        self.parameters['beta'] = s * ( 1 - mu )
#        
#        self.parameters['alpha'] = np.array( [
#                                              [1000, 10, 1],
#                                              [1000, 10, 1]
#                                              ], np.float )
#        
#        self.parameters['beta'] = np.array( [
#                                            [1, 10, 1000],
#                                            [1, 10, 1000]
#                                            ], np.float )
        
        print "Initial parameter values : ", self.parameters
    
    def _update_density_parameters( self ):        
        marginals = get_marginals( self.responsibilities, self.nclass )
        
        print "Begining numerical optimisation of alpha and beta."
        
        vars = []
        
        for genome in range( 2 ):
            for component in range( self.nclass ):
                a = self.data.a[genome]
                b = self.data.b[genome]
                
                x = np.zeros( ( 2, ) )
                
                x[0] = self.parameters['alpha'][genome][component]
                x[1] = self.parameters['beta'][genome][component]
                
                resp = marginals[genome][:, component]
                
                precision_prior = self.priors['precision'][genome][component]
                location_prior = self.priors['location'][genome][component]
                
                vars.append( [x, a, b, resp, precision_prior, location_prior] )
        
        results = []
        for var in vars:
            results.append( get_mle_p(var) )
        
#        results = self.pool.map( get_mle_p, vars )
        
        for genome in range( 2 ):
            for component in range( self.nclass ):
                i = genome * self.nclass + component
                
                self.parameters['alpha'][genome][component] = results[i][0]
                self.parameters['beta'][genome][component] = results[i][1]
                
class JointBinomialPosterior( EMPosterior ):
    def __init__( self, data, priors, responsibilities, nclass=3 ):
        self.nclass = nclass
        
        EMPosterior.__init__( self, data, priors, responsibilities )
    
    def _init_parameters( self ):
        self.parameters = {}
        
        self._update_mix_weights()
        
        self._update_density_parameters()
           
    def _update_density_parameters( self ):
        marginals = get_marginals( self.responsibilities, self.nclass )
        
        self.parameters['mu'] = []
        
        for genome in range( 2 ):
            a = self.data.a[genome]
            b = self.data.b[genome]
            
            alpha = self.priors['alpha'][genome]
            beta = self.priors['beta'][genome]
            
            tau = marginals[genome]
            
            self.parameters['mu'].append( self._update_mu( a, b, alpha, beta, tau ) )
    
    def _update_mu( self, a, b, alpha, beta, tau ):       
        d = a + b
        
        n = self.data.nrows
        shape = ( n, 1 )
        
        a = a.reshape( shape )
        d = d.reshape( shape )
        
        ref_sum = np.sum( tau * a, axis=0 )
        
        depth_sum = np.sum( tau * d, axis=0 )
        
        numerator = ref_sum + alpha - 1.
        
        denominator = depth_sum + alpha + beta - 2.
        
        return np.exp( np.log( numerator ) - np.log( denominator ) )
