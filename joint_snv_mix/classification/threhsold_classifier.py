'''
Created on 2011-02-10

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix import constants

class ThresholdModel( object ):
    def __init__( self, normal_threshold=0.05, tumour_threshold=0.1, min_var_depth=4 ):        
        self.threshold = {}
        self.threshold['normal'] = normal_threshold
        self.threshold['tumour'] = tumour_threshold

        self.min_var_depth = min_var_depth
    
    def classify( self, data ):
        genotypes = self._call_genotypes( data )
        
        joint_genotypes = self._call_joint_genotypes( data, genotypes )
        
        return joint_genotypes
    
    def _call_genotypes( self, data ):
        genotypes = {}

        n = data.a['normal'].size
        
        for genome in constants.genomes:
            a = data.a[genome]
            b = data.b[genome]
            d = np.asanyarray( a + b, dtype=np.float )
            
            freq = b / d
            t = self.threshold[genome]
            
            ref = ( freq < t )
            
            het = np.logical_and( t <= freq, freq <= ( 1 - t ) )
            
            hom = ( ( 1 - t ) < freq )
            
            genotypes[genome] = -1 * np.ones( ( n, ) )
            
            genotypes[genome][ref] = 0
            genotypes[genome][het] = 1
            genotypes[genome][hom] = 2
                    
        return genotypes
                    
    def _call_joint_genotypes( self, data, genotypes ):
        normal = genotypes['normal']
        tumour = genotypes['tumour']
        
        joint_genotypes = 3 * normal + tumour
        
        b_T = data.b['tumour']
        
        # Set low coverage sites to reference.
        joint_genotypes[b_T < self.min_var_depth] = 0

        return joint_genotypes
