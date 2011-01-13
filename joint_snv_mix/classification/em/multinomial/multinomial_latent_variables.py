'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2
from jsm_models.em.em_latent_variables import EMLatentVariables

class MultinomialLatentVariables( EMLatentVariables ):
    def _init_responsibilities( self, data ):
        '''
        Intialise responsibilities via k-means clustering.
        '''
        counts_1 = np.asarray( data.counts[0], dtype=np.float64 )
        counts_2 = np.asarray( data.counts[1], dtype=np.float64 )
        
        shape = ( counts_1.shape[0], 1 )
        
        p_1 = counts_1 / counts_1.sum( axis=1 ).reshape( shape )
        p_2 = counts_2 / counts_2.sum( axis=1 ).reshape( shape )
        
        shape = ( data.nrows, 100 )
        
        responsibilities = np.zeros( shape )
        
        init_centers = np.array( [
                                  ( 1., 0., 0., 0. ),
                                  ( 0.5, 0.5, 0., 0. ),
                                  ( 0.5, 0., 0.5, 0. ),
                                  ( 0.5, 0., 0., 0.5 ),
                                  ( 0., 1., 0., 0. ),
                                  ( 0., 0.5, 0.5, 0. ),
                                  ( 0., 0.5, 0., 0.5 ),
                                  ( 0., 0., 1., 0. ),
                                  ( 0., 0., 0.5, 0.5 ),
                                  ( 0., 0., 0., 1. )
                                  ] )
        
        cluster_centers_1, labels_1 = kmeans2( p_1, init_centers, minit='matrix' )
        cluster_centers_2, labels_2 = kmeans2( p_2, init_centers, minit='matrix' )

        labels = 10 * labels_1 + labels_2

        for id in range( 100 ):
            index = labels == id
            
            responsibilities[index, :] = 0.
            responsibilities[index, id] = 1.
        
        self.responsibilities = responsibilities / responsibilities.sum( axis=1 ).reshape( ( responsibilities.shape[0], 1 ) )
