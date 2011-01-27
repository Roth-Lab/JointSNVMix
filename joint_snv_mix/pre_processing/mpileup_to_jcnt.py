#!/usr/bin/env python
import tables
import warnings
warnings.filterwarnings( 'ignore', category=tables.NaturalNameWarning )

import csv

from collections import Counter

from joint_snv_mix.file_formats.pileup import parse_call_string

from joint_snv_mix.file_formats.jcnt import JointCountsFile

def main( args ):
    reader = get_reader( args.mpileup_file_name )
    jcnt_file = JointCountsFile( args.jcnt_file_name, 'w' )
    
    rows = {}
    i = 0
    
    for row in reader:
        chr_name = row['chr_name']
        
        if chr_name not in rows:
            rows[chr_name] = []
        
        chr_coord = int( row['chr_coord'] )
        
        normal_depth = int( row['normal_depth'] )
        tumour_depth = int( row['tumour_depth'] )
        
        ref_base = row['ref_base'].upper()
        
        # Skip lines below coverage threshold.
        if normal_depth < args.min_depth or tumour_depth < args.min_depth:
            continue
        
        normal_bases = parse_call_string( ref_base, row['normal_call_string'] )
        tumour_bases = parse_call_string( ref_base, row['tumour_call_string'] )
         
        normal_counter = Counter( normal_bases )
        tumour_counter = Counter( tumour_bases )
        
        normal_non_ref_base, normal_counts = parse_counts( normal_counter, ref_base )
        tumour_non_ref_base, tumour_counts = parse_counts( tumour_counter, ref_base )
        
        # Check for lines below reade depth. Necessary since * symbol counts against read depth above.
        d_N = normal_counts[0] + normal_counts[1]
        d_T = tumour_counts[0] + tumour_counts[1]
        if d_N < args.min_depth or d_T < args.min_depth:
            continue
        
        # Skip lines with no variants.
        if normal_counts[1] == 0 and tumour_counts[1] == 0:
            continue
        
        jcnt_entry = [ chr_coord, ref_base, normal_non_ref_base, tumour_non_ref_base ]
        jcnt_entry.extend( normal_counts )
        jcnt_entry.extend( tumour_counts )
        
        rows[chr_name].append( jcnt_entry )
        
        i += 1
        
        if i >= 1e4:
            print chr_name, chr_coord
            
            write_rows( jcnt_file, rows )
            
            rows = {}
            i = 0
    
    # Last call to write remaining rows.
    write_rows( jcnt_file, rows )
        
    jcnt_file.close()

def get_reader( mpileup_file_name ):
    csv.field_size_limit( 10000000 )
    
    fields = [
          'chr_name',
          'chr_coord',
          'ref_base',
          'normal_depth',
          'normal_call_string',
          'normal_base_qual_string',
          'tumour_depth',
          'tumour_call_string',
          'tumour_base_qual_string'
          ]

    reader = csv.DictReader( open( mpileup_file_name ), fieldnames=fields, delimiter='\t', quoting=csv.QUOTE_NONE )
        
    return reader

def parse_counts( counter, ref_base ):
    ref_counts = counter[ref_base]
    
    del counter[ref_base]
    
    # Check if there is any non-ref bases.
    if len( counter ) > 0:
        non_ref_base, non_ref_counts = counter.most_common( 1 )[0]
    else:
        non_ref_base = 'N'
        non_ref_counts = 0
    
    counts = ( ref_counts, non_ref_counts )
    
    return non_ref_base, counts

def write_rows( jcnt_file, rows ):
    for chr_name, chr_rows in rows.items():
        jcnt_file.add_rows( chr_name, chr_rows )
