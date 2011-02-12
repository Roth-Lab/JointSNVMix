'''
Created on 2011-02-11

@author: Andrew Roth
'''
from argparse import Namespace

from joint_snv_mix.file_formats.cnsm import ConanSnvMixReader
import csv
from joint_snv_mix import constants

def main( args ):
    reader = ConanSnvMixReader( args.cnsm_file_name )
    
    index_fields = [
                      'chrom',
                      'position',
                      'ref_base',
                      'normal_var_base',
                      'tumour_var_base',
                      'normal_counts_a',
                      'normal_counts_b',
                      'tumour_counts_a',
                      'tumour_counts_b'
                      ]
    
    cn_states = sorted( reader.get_cn_states() )
    
    for cn_state in  cn_states:
        chr_list = sorted( reader.get_chr_list( cn_state ) )
        
        out_file_name = args.prefix + "." + cn_state + ".tsv"
        
        writer = csv.writer( open( out_file_name, 'w' ), delimiter='\t' )
        
        som_indices = constants.conan_somatic_indices[cn_state]
        
        for chr_name in chr_list:
            index_rows, soft_labels = reader.get_rows( cn_state, chr_name )
            
            for i, row in enumerate( index_rows ):
                resp = soft_labels[i]
                                
                somatic_prob = resp[som_indices].sum()
                
                if somatic_prob >= 0.5:
                    print cn_state, chr_name, row[:], resp
                    row = row[:]
                    row = list( row )
                    row.insert( 0, chr_name )
                    row.extend( resp )
                    
                    writer.writerow( row )
    
    reader.close()         
                    
if __name__ == "__main__":
    import sys
    
    args = Namespace()
    
    args.cnsm_file_name = sys.argv[1]
    args.prefix = sys.argv[2]
    
    main( args )
