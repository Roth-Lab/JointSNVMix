'''
Created on 2011-01-18

@author: Andrew Roth
'''
import numpy as np
import multiprocessing
from joint_snv_mix.classification.utils.beta_binomial_map_estimators import get_mle_p, get_junk_map
from joint_snv_mix import constants

def get_marginals( responsibilities, nclass ):
    nrows = responsibilities.shape[0]

    shape = ( nrows, nclass, nclass )
    
    responsibilities = responsibilities.reshape( shape )
    
    marginals = {}
    
    marginals['normal'] = responsibilities.sum( axis=2 )
    marginals['tumour'] = responsibilities.sum( axis=1 )

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
            
            precision_prior = self.priors['precision']
            location_prior = self.priors['location']
            
            vars.append( [x, a, b, resp, location_prior, precision_prior, component] )
                
#        results = self.pool.map( get_mle_p, vars )
        
        results = []        
        for var in vars:
            results.append( get_mle_p( var ) )
        
        for component in range( 3 ):
            self.parameters['alpha'][component] = results[component][0]
            self.parameters['beta'][component] = results[component][1]

class IndependentBinomialPosterior( EMPosterior ):
    def _init_parameters( self ):
        self.parameters = {}
        
        self._update_mix_weights()

        alpha = self.priors['mu']['alpha']
        beta = self.priors['mu']['beta']
        
        self.parameters['mu'] = alpha / ( alpha + beta )
        
        print self.parameters

    def _update_density_parameters( self ):
        a = self.data.a
        b = self.data.b
        d = a + b
        
        shape = ( self.data.nrows, 1 )
        a = a.reshape( shape )
        d = d.reshape( shape )
        
        alpha = self.priors['mu']['alpha']
        beta = self.priors['mu']['beta']
        
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
               
        EMPosterior.__init__( self, data, priors, responsibilities )
    
    def _init_parameters( self ):
        '''
        Initialise parameters. This is only necessary to initialise gradient descent. 
        '''
        self.parameters = {}

        self._update_mix_weights()

        for genome in constants.genomes:
            self.parameters[genome] = {}
            
            location_alpha = self.priors[genome]['location']['alpha']
            location_beta = self.priors[genome]['location']['beta']
        
            precision_shape = self.priors[genome]['precision']['shape']
            precision_scale = self.priors[genome]['precision']['scale'] 
        
            s = precision_shape * precision_scale
            mu = location_alpha / ( location_alpha + location_beta ) 
        
            self.parameters[genome]['alpha'] = s * mu
            self.parameters[genome]['beta'] = s * ( 1 - mu )
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
        
        for genome in constants.genomes:
            for component in range( self.nclass ):
                a = self.data.a[genome]
                b = self.data.b[genome]
                
                x = np.zeros( ( 2, ) )
                
                x[0] = self.parameters[genome]['alpha'][component]
                x[1] = self.parameters[genome]['beta'][component]
                
                resp = marginals[genome][:, component]
                
                precision_prior = self.priors[genome]['precision']
                location_prior = self.priors[genome]['location']
                
                vars.append( [x, a, b, resp, location_prior, precision_prior, component] )        
        
#        results = []
#        for var in vars:
#            results.append( get_mle_p( var ) )

        ncpus = self.nclass * 2
        
        pool = multiprocessing.Pool( processes=ncpus, maxtasksperchild=1 )
        results = pool.map( get_mle_p, vars )
        pool.close()
                
        for i, genome in enumerate( constants.genomes ):
            for component in range( self.nclass ):
                index = i * self.nclass + component
                
                self.parameters[genome]['alpha'][component] = results[index][0]
                self.parameters[genome]['beta'][component] = results[index][1]
                
        self.parameters['junk'] = get_junk_map( self.data, self.parameters['junk'], self.responsibilities[:, 10] )
        
                
class JointBinomialPosterior( EMPosterior ):
    def __init__( self, data, priors, responsibilities, nclass=3 ):
        self.nclass = nclass
        
        EMPosterior.__init__( self, data, priors, responsibilities )
    
    def _init_parameters( self ):
        self.parameters = {}
        self.parameters['normal'] = {}
        self.parameters['tumour'] = {}
        
        self._update_mix_weights()
        
        self._update_density_parameters()
           
    def _update_density_parameters( self ):
        marginals = get_marginals( self.responsibilities, self.nclass )
        
        for genome in constants.genomes:
            a = self.data.a[genome]
            b = self.data.b[genome]
            
            alpha = self.priors[genome]['mu']['alpha']
            beta = self.priors[genome]['mu']['beta']
            
            tau = marginals[genome]
            
            self.parameters[genome]['mu'] = self._update_mu( a, b, alpha, beta, tau )
    
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
