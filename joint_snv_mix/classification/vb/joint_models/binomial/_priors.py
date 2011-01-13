'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

class _Priors( object ):
    def __init__( self ):
#        self.kappa = np.array([1e5, 10, 5, 10, 1000, 10, 1, 1, 100])
        self.kappa = np.array( [1.] * 9 ) * 10
        
        self.alpha_1 = np.array( [99., 5., 0.2] )
        self.alpha_2 = np.array( [99., 5., 0.2] )

        self.beta_1 = np.array( [0.2, 5., 99.] )
        self.beta_2 = np.array( [0.2, 5., 99.] )
