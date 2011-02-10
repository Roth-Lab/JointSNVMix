'''
Created on 2011-02-09

@author: Andrew Roth
'''
import numpy as np

from fisher import pvalue_npy

from joint_snv_mix import constants

class FisherModel( object ):
    def __init__( self, p_value_threshold=0.05, base_line_error=0.001,
                  min_var_freq=0.1, min_hom_freq=0.8, min_var_support=2 ):
        
        self.p_value_threshold = p_value_threshold
        
        self.base_line_error = base_line_error
        
        self.min_var_freq = min_var_freq
        self.min_hom_freq = min_hom_freq
        self.min_var_depth = min_var_support
    
    def classify( self, data ):
        genotypes = self._call_genotypes( data )
        
        joint_genotypes, germline_scores, somatic_scores = self._call_joint_genotypes( data, genotypes )
        
        return joint_genotypes, germline_scores, somatic_scores
    
    def _call_genotypes( self, data ):
        genotypes = {}

        n = data.a['normal'].size
        
        for genome in constants.genomes:
            a = data.a[genome]
            b = data.b[genome]
            d = np.asanyarray( a + b, dtype=np.float )
                        
            p_values = self._get_significance( a, b, self.base_line_error )
                        
            below_threshold = ( p_values <= self.p_value_threshold )
            
            var_freq = b / d            

            # Find variants below p-value threshold and above minimum variant frequency.
            above_min_var_freq = ( var_freq >= self.min_var_freq )
            
            var_indices = np.logical_and( below_threshold, above_min_var_freq )
            
            # Call homozygous if above minimum threshold.
            above_min_hom_freq = ( var_freq >= self.min_hom_freq )                        
            hom_var_indices = np.logical_and( var_indices, above_min_hom_freq )
            
            not_above_min_hom_freq = np.logical_not( above_min_hom_freq )
            het_var_indices = np.logical_and( var_indices, not_above_min_hom_freq )
            
            # Apply min variant depth filter
            above_min_var_depth = ( b >= self.min_var_depth )
            
            het_var_indices = np.logical_and( het_var_indices, above_min_var_depth )
            hom_var_indices = np.logical_and( hom_var_indices, above_min_var_depth )
            
            # Assign genotypes assuming reference as default
            genotypes[genome] = np.zeros( ( n, ) )
            
            genotypes[genome][het_var_indices] = 1
            genotypes[genome][hom_var_indices] = 2
                    
        return genotypes
                    
    def _call_joint_genotypes( self, data, genotypes ):
        normal = genotypes['normal']
        tumour = genotypes['tumour']
        
        normal_aa = ( normal == 0 )
        normal_ab = ( normal == 1 )
        normal_bb = ( normal == 2 )
            
        normal_var = np.logical_or( normal_ab, normal_bb )
        
        tumour_aa = ( tumour == 0 )
        tumour_ab = ( tumour == 1 )
        tumour_bb = ( tumour == 2 )
        
        tumour_var = np.logical_or( tumour_ab, tumour_bb )
                
        reference = np.logical_and( normal_aa, tumour_aa )
        germline = np.logical_and( normal_var, tumour_var )
        
        a_N = data.a['normal']
        b_N = data.b['normal']
        d_N = np.asanyarray( a_N + b_N, dtype=np.float )      
        normal_freq = b_N / d_N
        
        a_T = data.a['tumour']
        b_T = data.b['tumour']
        
        p_values = self._get_significance( a_T, b_T, normal_freq )
        
        non_match = ( normal != tumour )
        significant_p_values = ( p_values <= self.p_value_threshold )        
        significant_non_match = np.logical_and( non_match, significant_p_values )
        
        somatic = np.logical_and( normal_aa, significant_non_match )
        loh = np.logical_and( normal_ab, significant_non_match )
        
        uknown = np.logical_and( normal_bb, significant_non_match )
        
        non_significant_p_values = np.logical_not( significant_p_values )
        non_significant_non_match = np.logical_and( non_match, non_significant_p_values )        
        germline = np.logical_or( germline, non_significant_non_match )
        
        a = a_N + a_T
        b = b_N + b_T
        
        germline_pvalues = self._get_significance( a, b, self.base_line_error )
        
        n = a_N.size
        joint_genotypes = -1 * np.ones( ( n, ) )
        
        joint_genotypes[reference] = 0
        joint_genotypes[germline] = 1
        joint_genotypes[somatic] = 2
        joint_genotypes[loh] = 3
        joint_genotypes[uknown] = 4
        
        return joint_genotypes, germline_pvalues, p_values

    def _get_significance( self, a, b, expected_freq ):
        '''
        Compute the p-value where the null hypothesis is that b was obtained due to error.
        
        a, b : numpy array of counts
        base_line_error : expected error rate.
        '''
        d = a + b
        
        expected_b = np.around( d * expected_freq )
        
        expected_a = d - expected_b
        
        # Downcast to uint to work with fisher exact test function.
        a = np.asarray( a, dtype=np.uint )
        b = np.asarray( b, dtype=np.uint )
        expected_a = np.asarray( expected_a, dtype=np.uint )
        expected_b = np.asarray( expected_b, dtype=np.uint )
        
        left_tail, right_tail, two_tail = pvalue_npy( expected_a,
                                                      expected_b,
                                                      a,
                                                      b )
        
        return right_tail
