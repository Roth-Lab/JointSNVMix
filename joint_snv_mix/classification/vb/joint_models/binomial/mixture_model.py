'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from ._ess import _ESS
from ._posterior import _Posterior
from ._latent_variables import _LatentVariables
from ._lower_bound import _LowerBound

class MixtureModel:
    def train( self, data, priors, max_iter, threshold ):
        lv = _LatentVariables( data )
        ess = _ESS( data )
        posterior = _Posterior( priors )
        
        lower_bound = _LowerBound( priors, posterior, ess, lv )
        
        lb_change = float( 'inf' )
        old_lb = float( '-inf' )

        iter = 0

        while lb_change > threshold and iter < max_iter:
            lv.update( posterior )
            ess.update( lv )
            posterior.update( ess )

            lb = lower_bound.get_lower_bound()

            lb_change = lb - old_lb
            old_lb = lb

            if lb_change < 0:
                print "Warning posterior decreased. Exiting EM optimisation and classifying data points."
                print "Posterior change ", lb_change                    
                print posterior.kappa_bar / posterior.kappa_bar.sum()
                print posterior.alpha_bar_1 / ( posterior.alpha_bar_1 + posterior.beta_bar_1 )
                print posterior.alpha_bar_2 / ( posterior.alpha_bar_2 + posterior.beta_bar_2 )
                break

            print posterior.kappa_bar / posterior.kappa_bar.sum()
            print posterior.alpha_bar_1 / ( posterior.alpha_bar_1 + posterior.beta_bar_1 )
            print posterior.alpha_bar_2 / ( posterior.alpha_bar_2 + posterior.beta_bar_2 )

            iter += 1
            print iter, lb
            print lb_change
            
        parameters = {}
        parameters['pi'] = posterior.kappa_bar / posterior.kappa_bar.sum()
        mu_1 = posterior.alpha_bar_1 / ( posterior.alpha_bar_1 + posterior.beta_bar_1 )
        mu_2 = posterior.alpha_bar_2 / ( posterior.alpha_bar_2 + posterior.beta_bar_2 )
        parameters['mu'] = np.vstack( ( mu_1, mu_2 ) )

        return posterior
        

    def classify( self, data, posterior ):
        lv = _LatentVariables( data )
        
        lv.update( posterior )

        return lv.responsibilities
