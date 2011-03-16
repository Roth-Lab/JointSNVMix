#! /usr/bin/env python
'''
Created on 2010-08-10

@author: Andrew Roth
'''
import csv

import joint_snv_mix.constants as constants

from joint_snv_mix.file_formats.jsm import JointSnvMixReader

def main( args ):
    reader = JointSnvMixReader( args.jsm_file_name )
    writer = csv.writer( open( args.call_file_name, 'w' ), delimiter="\t" )

    header = get_header()
    writer.writerow( header )

    chr_list = sorted( reader.get_table_list() )

    for chr_name in sorted( chr_list ):
        if args.argmax:
            rows = reader.get_genotype_rows_by_argmax( chr_name, args.genotype_class )
        else:
            rows = reader.get_genotype_rows_by_prob( chr_name, args.genotype_class, args.prob_threshold )
        
        if not len(rows) == 0:
            continue
        
        rows = rows.tolist()
        
        rows = [list( row ) for row in rows]
        
        for row in rows:
            row[0] = "{0}:{1}".format( chr_name, row[0] )
        
        writer.writerows( rows )

    reader.close()

def get_header():
    header = []
    header.append( 'position' )
    header.append( 'ref_base' )
    header.append( 'normal_base' )
    header.append( 'tumour_base' )
    header.append( 'normal_ref_counts' )
    header.append( 'normal_non_ref_counts' )
    header.append( 'tumour_ref_counts' )
    header.append( 'tumour_non_ref_counts' )

    for genotype in constants.joint_genotypes:
        header.append( "p" + str( genotype ) )
    
    return header

