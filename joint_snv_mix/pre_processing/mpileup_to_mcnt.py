#!/usr/bin/env python
import tables
import warnings
from joint_snv_mix.constants import nucleotides
import bz2
warnings.filterwarnings( 'ignore', category=tables.NaturalNameWarning )

import csv

from collections import Counter

from joint_snv_mix.file_formats.pileup import parse_call_string

from joint_snv_mix.file_formats.mcnt import MultinomialCountsFile

def main( args ):
    if args.bzip2:
        mpileup_file = bz2.BZ2File( args.mpileup_file_name )
    else:
        mpileup_file = open( args.mpileup_file_name )
    
    reader = get_reader( mpileup_file )
    mcnt_file = MultinomialCountsFile( args.mcnt_file_name, 'w' )
    
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
        
        normal_counts = []
        tumour_counts = []
    
        variant = True
    
        for nucleotide in nucleotides:
            nc = normal_counter[nucleotide]
            tc = tumour_counter[nucleotide]
            
            if nucleotide != ref_base:
                if nc > 0 or tc > 0:
                    variant = True
            
            normal_counts.append( nc )
            tumour_counts.append( tc )
        
        if sum( normal_counts ) < args.min_depth or sum( tumour_counts ) < args.min_depth:
            continue
        
        # Skip lines with no variants.
        if not variant:
            continue
        
        mcnt_entry = [ chr_coord, ref_base ]
        mcnt_entry.extend( normal_counts )
        mcnt_entry.extend( tumour_counts )
        
        rows[chr_name].append( mcnt_entry )
        
        i += 1
        
        if i >= 1e4:
            print chr_name, chr_coord, mcnt_entry
            
            write_rows( mcnt_file, rows )
            
            rows = {}
            i = 0
    
    # Last call to write remaining rows.
    write_rows( mcnt_file, rows )
        
    mcnt_file.close()

def get_reader( mpileup_file ):
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

    reader = csv.DictReader( mpileup_file, fieldnames=fields, delimiter='\t', quoting=csv.QUOTE_NONE )
        
    return reader

def write_rows( mcnt_file, rows ):
    for chr_name, chr_rows in rows.items():
        mcnt_file.add_rows( chr_name, chr_rows )
