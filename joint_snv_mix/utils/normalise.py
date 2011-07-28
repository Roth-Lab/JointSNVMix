'''
Created on 2010-11-22

@author: Andrew Roth
'''
import numpy as np

def log_space_normalise_rows( log_X ):
    nrows = log_X.shape[0]
    shape = ( nrows, 1 )
    
    log_norm_const = np.logaddexp.reduce( log_X, axis=1 )
    log_norm_const = log_norm_const.reshape( shape )

    log_X = log_X - log_norm_const
    
    X = np.exp( log_X )
    
    dt = X.dtype
    eps = np.finfo( dt ).eps
    
    X[X <= eps] = 0.
    
    return X
