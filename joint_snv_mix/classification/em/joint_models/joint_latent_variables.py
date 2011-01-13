'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2

from joint_snv_mix.classification.em.em_latent_variables import EMLatentVariables


class JointLatentVariables( EMLatentVariables ):
    def _init_responsibilities( self, data ):
        '''
        Intialise responsibilities via k-means clustering.
        '''
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
        
        self.responsibilities = responsibilities
