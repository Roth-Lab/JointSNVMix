'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2

from jsm_models.em.em_latent_variables import EMLatentVariables

class IndependentLatenVariables( EMLatentVariables ):        
    def _init_responsibilities( self, data ):
        a = np.asarray( data.a, dtype=np.float64 )
        b = np.asarray( data.b, dtype=np.float64 )
        
        shape = ( data.nrows, 3 )
        
        responsibilities = np.zeros( shape )
        
        p = a / ( a + b )
        
        init_centers = np.array( [0.99, 0.5, 0.01] )
        
        cluster_centers, labels = kmeans2( p, init_centers )
        
        sorted_centers = np.argsort( cluster_centers )
        
        for id in sorted_centers:
            index = labels == id
            
            responsibilities[index, id] = 1.0
        
        self.responsibilities = responsibilities
