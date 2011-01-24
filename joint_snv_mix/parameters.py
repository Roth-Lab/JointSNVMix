'''
Created on 2011-01-21

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix.classification.utils.log_pdf import log_gamma_pdf, log_beta_pdf
from joint_snv_mix.classification.utils.beta_binomial_map_estimators import get_mle_p
import multiprocessing

class Parameters( object ):
    pass

class BetaBinomialParameters( Parameters ):
    def __init__( self, precision_shape, precision_scale, location_alpha, location_beta ):
        self.precision_shape = precision_shape 
        self.precision_scale = precision_scale
        
        self.location_alpha = location_alpha
        self.location_beta = location_beta
    
        self.nclass = self.location_alpha.size
    
        self.pool = multiprocessing.Pool( processes=self.nclass, maxtasksperchild=1 )
        
        self._init_values()
        
    def _init_values( self ):
        '''
        Initialise parameters. This is only necessary to initialise gradient descent. 
        '''        
        values = np.linspace( 1, 99, self.nclass )
   
        self.alpha = np.vstack( ( values[::-1], values[::-1] ) )
        
        self.beta = np.vstack( ( values, values ) )
        
        print "Initial parameter values : ", self.alpha, self.beta
        
    def get_pdf( self ):
        s = self.alpha + self.beta
        mu = self.alpha / s
    
        precision_term = log_gamma_pdf( s, self.precision_shape, self.precision_scale )
        
        location_term = log_beta_pdf( mu, self.location_alpha, self.location_beta )
                
        return precision_term + location_term
    
    def update( self, a, b, responsibilities ):
        print "Begining numerical optimisation of alpha and beta."
        
        nclass = responsibilities.shape[1]
        
        vars = []
    
        for component in range( nclass ):           
            x = np.zeros( ( 2, ) )
            
            x[0] = self.alpha[component]
            x[1] = self.beta[component]
            
            resp = responsibilities[:, component]
            
            precision_prior = [ self.precision_shape[component],
                                self.precision_scale[component] ]
            
            location_prior = [ self.location_alpha[component],
                               self.location_alpha[component] ]
            
            vars.append( [x, a, b, resp, component, precision_prior, location_prior] )
    
        results = self.pool.map( get_mle_p, vars )
        
        
        for component in range( self.nclass ):        
            self.alpha[component] = results[component][0]
            self.beta[component] = results[component][1]
