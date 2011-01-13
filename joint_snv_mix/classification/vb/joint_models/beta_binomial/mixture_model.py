'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2

from ._posterior import _Posterior
from ._latent_variables import _LatentVariables
from ._lower_bound import _LowerBound
from ...._data import _JointData
from models.utils.beta_binomial_parameter_estimators import get_mle_p
import multiprocessing

class MixtureModel:
    def __init__( self ):
        self.priors = _Priors()

    def train( self, X, threshold=1e-8, max_iter=100 ):
        data = _JointData( X )
        
        lv = _LatentVariables( data, self.priors )
        posterior = _Posterior( self.priors )
        
        lv.set_posterior( posterior )
        posterior.set_latent_variables( lv )
                
        lower_bound = _LowerBound( self.priors, lv, posterior )
        
        lb_change = float( 'inf' )
        old_lb = float( '-inf' )

        iter = 0
        
        repsonsibilities = self._init_responsibilities( data ) 
        lv.set_responsibilities( repsonsibilities )

        while abs( lb_change ) > threshold and iter < max_iter:            
            posterior.update()
#            self.update_alpha_beta( data, lv, posterior )
            
            lv.update()
            
            lb = lower_bound.get_lower_bound()

            lb_change = lb - old_lb
            old_lb = lb

#            if lb_change < 0:
#                print "Warning posterior decreased. Exiting EM optimisation and classifying data points."
#                print "Posterior change ", lb_change                    
#
#                print posterior.kappa_bar / posterior.kappa_bar.sum()
#                print posterior.parameters.alpha
#                print posterior.parameters.beta
#
#                break

            print posterior.kappa_bar / posterior.kappa_bar.sum()
            print posterior.parameters.alpha
            print posterior.parameters.beta

            iter += 1
            print iter, lb
            print lb_change

        parameters = {}
        parameters['pi'] = posterior.kappa_bar / posterior.kappa_bar.sum()
        
        return parameters, lv.responsibilities
    
    def update_alpha_beta( self, data, lv, posterior ):
        
        vars = []
        for genome in range( 2 ):
            a = data.a[genome]
            b = data.b[genome]
            
            params = zip( posterior.parameters.alpha[genome], posterior.parameters.beta[genome] )
            
            for component, x in enumerate( params ):
                resp = lv.marginals[genome][:, component]
                
                vars.append( [x, a, b, resp] )
        
        p = multiprocessing.Pool()
        
        results = p.map( get_mle_p, vars )
        
        for genome in range( 2 ):
            for component in range( 3 ):
                i = genome * 3 + component
                
                posterior.parameters.alpha[genome][component] = results[i][0]
                posterior.parameters.beta[genome][component] = results[i][1]
            
        
    def _init_responsibilities( self, data ):
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
            
            responsibilities[index, :] = 0.
            responsibilities[index, id] = 1.
        
        return responsibilities
    
class _Priors( object ):
    def __init__( self ):
        self.kappa = np.array( [100.] * 9 )
