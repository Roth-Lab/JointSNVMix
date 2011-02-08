'''
Created on 2010-09-15

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2

from joint_snv_mix.classification.utils.normalise import log_space_normalise_rows
from joint_snv_mix.classification.likelihoods import independent_beta_binomial_log_likelihood, \
    independent_binomial_log_likelihood, joint_beta_binomial_log_likelihood, joint_binomial_log_likelihood, \
    joint_multinomial_log_likelihood

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

class IndependentLatenVariables( EMLatentVariables ):        
    def _init_responsibilities( self, data ):
        a = np.asarray( data.a, dtype=np.float64 )
        b = np.asarray( data.b, dtype=np.float64 )
        
        shape = ( data.nrows, 3 )
        
        responsibilities = np.zeros( shape )
        
        p = a / ( a + b )
        
        init_centers = np.array( [1., 0.5, 0.] )
        
        cluster_centers, labels = kmeans2( p, init_centers )
        
        sorted_centers = np.argsort( cluster_centers )
        
        for id in sorted_centers:
            index = labels == id
            
            responsibilities[index, id] = 1.0
        
        self.responsibilities = responsibilities

class JointLatentVariables( EMLatentVariables ):
    def _init_responsibilities( self, data ):
        '''
        Intialise responsibilities via k-means clustering.
        '''
        a_1 = np.asarray( data.a['normal'], dtype=np.float64 )
        b_1 = np.asarray( data.b['normal'], dtype=np.float64 )
        p_1 = a_1 / ( a_1 + b_1 )
              
        a_2 = np.asarray( data.a['tumour'], dtype=np.float64 )
        b_2 = np.asarray( data.b['tumour'], dtype=np.float64 )
        p_2 = a_2 / ( a_2 + b_2 )

        shape = ( data.nrows, 9 )
        
        responsibilities = np.zeros( shape )
        
        init_centers = np.array( ( 1., 0.5, 0. ) )
        
        cluster_centers_1, labels_1 = kmeans2( p_1, init_centers, minit='matrix' )
        cluster_centers_2, labels_2 = kmeans2( p_2, init_centers, minit='matrix' )

        labels = 3 * labels_1 + labels_2

        for id in range( 9 ):
            index = labels == id
            
            responsibilities[index, id] = 1.
        
        self.responsibilities = responsibilities

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
        
        self.responsibilities = responsibilities

#=======================================================================================================================
# Independent Models
#=======================================================================================================================
class IndependentBetaBinomialLatentVariables( IndependentLatenVariables ):
    def __init__( self, data ):
        IndependentLatenVariables.__init__( self, data )
        
        self.likelihood_func = independent_beta_binomial_log_likelihood
        
class IndependentBinomialLatentVariables( IndependentLatenVariables ):
    def __init__( self, data ):
        IndependentLatenVariables.__init__( self, data )
        
        self.likelihood_func = independent_binomial_log_likelihood

#=======================================================================================================================
# Joint Models
#=======================================================================================================================
class JointBetaBinomialLatentVariables( JointLatentVariables ):
    def __init__( self, data ):
        JointLatentVariables.__init__( self, data )
        
        self.likelihood_func = joint_beta_binomial_log_likelihood
        
class JointBinomialLatentVariables( JointLatentVariables ):
    def __init__( self, data ):
        JointLatentVariables.__init__( self, data )
        
        self.likelihood_func = joint_binomial_log_likelihood
        
#=======================================================================================================================
# Multinomial
#=======================================================================================================================
class JointMultinomialLatentVariables( MultinomialLatentVariables ):
    def __init__( self, data ):
        MultinomialLatentVariables.__init__( self, data )
        
        self.likelihood_func = joint_multinomial_log_likelihood
