'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

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
