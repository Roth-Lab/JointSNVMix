'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2
from jsm_models.em.em_latent_variables import EMLatentVariables

nclass = 15

class ExtendedMultinomialLatentVariables( EMLatentVariables ):
    def _init_responsibilities( self, data ):
        '''
        Intialise responsibilities via k-means clustering.
        '''
        counts_1 = np.asarray( data.counts[0], dtype=np.float64 )
        counts_2 = np.asarray( data.counts[1], dtype=np.float64 )
        
        shape = ( counts_1.shape[0], 1 )
        
        p_1 = counts_1 / counts_1.sum( axis=1 ).reshape( shape )
        p_2 = counts_2 / counts_2.sum( axis=1 ).reshape( shape )
        
        shape = ( data.nrows, nclass ** 2 )
        
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
                                  ( 0.33, 0.33, 0.33, 0. ),
                                  ( 0.33, 0., 0.33, 0.33 ),
                                  ( 0., 0.33, 0.33, 0.33 ),
                                  ( 0.25, 0.25, 0.25, 0.25 )
                                  ] )
        
        cluster_centers_1, labels_1 = kmeans2( p_1, init_centers, minit='matrix' )
        cluster_centers_2, labels_2 = kmeans2( p_2, init_centers, minit='matrix' )

        labels = nclass * labels_1 + labels_2

        for id in range( nclass ** 2 ):
            index = labels == id
            
            responsibilities[index, :] = 0.
            responsibilities[index, id] = 1.
        
        self.responsibilities = responsibilities / responsibilities.sum( axis=1 ).reshape( ( responsibilities.shape[0], 1 ) )
