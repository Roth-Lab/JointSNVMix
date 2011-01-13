'''
Created on 2010-12-11

Port of Varscan.class from Varscan 2.2. 

@author: Andrew Roth
'''
import csv

import numpy as np

from fisher import pvalue_npy
from joint_snv_mix.file_formats.jcnt import JointCountsReader
from jsm_models.data import JointData

def main( jcnt_file_name, tsv_file_name ):
    reader = JointCountsReader( jcnt_file_name )
    writer = csv.writer( open( tsv_file_name, 'w' ), delimiter='\t' )
    
    chr_list = reader.get_chr_list()
    
    for chr_name in sorted( chr_list ):
        counts = reader.get_counts( chr_name )
        
        data = JointData( counts )
        
        labels, scores = classify( data.a[0], data.b[0], data.a[1], data.b[1], 0.05 )
              
        jcnt_rows = reader.get_rows( chr_name )
        
        write_rows( writer, chr_name, labels, scores, jcnt_rows )
    
    reader.close()

def write_rows( writer, chr_name, labels, scores, jcnt_rows ):
    rows = []
    for jcnt_row, score, label in zip( jcnt_rows, scores, labels ):
        row = [chr_name]
        row.extend( jcnt_row.tolist() )
        row.extend( score.tolist() )
        row.append( label )
        
        rows.append( row )
    
    writer.writerows( rows )

def classify( normal_a, normal_b, tumour_a, tumour_b, p_value_threshold ):
    # Check for variant positions in tumour
    normal_is_variant, normal_p_values = call_variants( normal_a, normal_b, p_value_threshold )
    tumour_is_variant, tumour_p_values = call_variants( tumour_a, tumour_b, p_value_threshold )
    
    normal_is_not_variant = np.logical_not( normal_is_variant )
    tumour_is_not_variant = np.logical_not( tumour_is_variant )

    reference_indices = np.logical_and( normal_is_not_variant, tumour_is_not_variant )
    germline_indices = np.logical_and( normal_is_variant, tumour_is_variant )
    somatic_indices = np.logical_and( normal_is_not_variant, tumour_is_variant )
    loh_indices = np.logical_and( normal_is_variant, tumour_is_not_variant )
    
    n = normal_a.size
    labels = np.zeros( ( n, ) , dtype=np.uint8 )
    
    # Reference = 0, Germline = 1, Somatic = 2, LOH = 3 
    labels[reference_indices] = 0
    labels[germline_indices] = 1
    labels[somatic_indices] = 2
    labels[loh_indices] = 3
    
    normal_variant_score = 1 - normal_p_values
    tumour_variant_score = 1 - tumour_p_values
    
    scores = np.zeros( ( n, 4 ), dtype=np.float64 )
    
    scores[:, 0] = ( 1 - normal_variant_score ) * ( 1 - tumour_variant_score )
    scores[:, 1] = normal_variant_score * tumour_variant_score
    scores[:, 2] = ( 1 - normal_variant_score ) * tumour_variant_score
    scores[:, 3] = normal_variant_score * ( 1 - tumour_variant_score )
    
    return labels, scores

def call_variants( a, b, max_p_value ):
    '''
    Return numpy array of booleans indicating if position is a variant.
    '''
    p_values = get_significance_vs_baseline( a, b )
    
    is_variants = p_values <= max_p_value
    
    return is_variants, p_values

def get_significance_vs_baseline( a, b, base_line_error=0.001 ):
    '''
    Compute the p-value where the null hypothesis is that b was obtained due to error.
    
    a, b : numpy array of counts
    base_line_error : expected error rate.
    '''
    d = a + b
    
    expected_b = np.around( d * base_line_error )
    
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

if __name__ == "__main__":
    import sys
    
    jcnt_file_name = sys.argv[1]
    tsv_file_name = sys.argv[2]

    main( jcnt_file_name, tsv_file_name )
