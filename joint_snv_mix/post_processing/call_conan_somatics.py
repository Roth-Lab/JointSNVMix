'''
Created on 2011-02-11

@author: Andrew Roth
'''
from argparse import Namespace

from joint_snv_mix.file_formats.cnsm import ConanSnvMixReader
import csv
from joint_snv_mix import constants

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

def main( args ):
    reader = ConanSnvMixReader( args.cnsm_file_name )
    

    
    cn_states = sorted( reader.get_cn_states() )
    
    for cn_state in  cn_states:
        chr_list = sorted( reader.get_chr_list( cn_state ) )
        
        writer, p_genotype_str = get_writer( args, cn_state )
        
        som_indices = constants.conan_somatic_indices[cn_state]
        
        for chr_name in chr_list:
            index_rows, soft_labels = reader.get_rows( cn_state, chr_name )
            
            for i, index_row in enumerate( index_rows ):
                resp = soft_labels[i]
                                
                somatic_prob = resp[som_indices].sum()
                
                if somatic_prob >= 0.5:
                    row = format_row( chr_name, index_row, resp, p_genotype_str )
                    
                    writer.writerow( row )
    
    reader.close()
 
def get_writer( args, cn_state ):
    out_file_name = args.prefix + "." + cn_state + ".tsv"
    
    fields = index_fields[:]
    
    p_genotype_str = [
                        "_".join( ( 'p', x[0], x[1] ) ) 
                        for x in constants.conan_joint_genotypes[cn_state]
                        ]
    
    fields.extend( p_genotype_str )
    
    writer = csv.DictWriter( open( out_file_name, 'w' ), fields, delimiter='\t' )
    
    writer.writeheader()
    
    return writer, p_genotype_str
    
def format_row( chr_name, index_row, resp, p_genotype_str ):
    row = {}
    
    index_row = list( index_row[:] )
    
    for i, field_name in enumerate( index_fields ):
        row[field_name] = index_row[i - 1]
    
    for p_genotype in p_genotype_str:
        row[field_name] = row[p_genotype]
        
    return row
                    
if __name__ == "__main__":
    import sys
    
    args = Namespace()
    
    args.cnsm_file_name = sys.argv[1]
    args.prefix = sys.argv[2]
    
    main( args )
