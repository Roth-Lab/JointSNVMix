'''
Created on 2010-09-15

@author: Andrew Roth
'''
from jsm_models.utils.normalise import log_space_normalise_rows

class EMLatentVariables( object ):
    def __init__( self, data ):       
        self.data = data
        
        self._init_responsibilities( data )
        
        self.likelihood_func = None

    def update( self, parameters ):       
        log_responsibilities = self.likelihood_func( self.data, parameters )

        self.responsibilities = log_space_normalise_rows( log_responsibilities )
       
    def _init_responsibilities( self, data ):
        NotImplementedError
