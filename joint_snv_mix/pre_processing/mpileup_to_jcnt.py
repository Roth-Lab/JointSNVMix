#!/usr/bin/env python
import bz2
import csv

from collections import Counter

import numpy as np

import tables
import warnings
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

from joint_snv_mix.file_formats.pileup import parse_call_string

from joint_snv_mix.file_formats.jcnt import JointCountsFile

ascii_offset = 33

def main(args):
    if args.bzip2:
        mpileup_file = bz2.BZ2File(args.mpileup_file_name)
    else:
        mpileup_file = open(args.mpileup_file_name)
    
    reader = get_reader(mpileup_file)
    
    jcnt_file = JointCountsFile(args.jcnt_file_name, 'w')
    
    rows = {}
    i = 0
    
    for row in reader:
        chr_name = row['chr_name']
    
        if chr_name not in rows:
            rows[chr_name] = []
        
        jcnt_entry = get_jcnt_entry(row, args.min_depth, args.min_qual)
        
        if jcnt_entry is None:
            continue
        
        rows[chr_name].append(jcnt_entry)
        
        i += 1
        
        if i >= 1e5:
            print rows[-1]
            
            write_rows(jcnt_file, rows)
            
            rows = {}
            i = 0
    
    # Last call to write remaining rows.
    write_rows(jcnt_file, rows)
        
    jcnt_file.close()
    mpileup_file.close()

def get_reader(mpileup_file):
    csv.field_size_limit(10000000)
    
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

    reader = csv.DictReader(mpileup_file, fieldnames=fields, delimiter='\t', quoting=csv.QUOTE_NONE)
        
    return reader

def get_jcnt_entry(row, min_depth, min_qual):    
    chr_coord = int(row['chr_coord'])
    
    normal_depth = int(row['normal_depth'])
    tumour_depth = int(row['tumour_depth'])
    
    ref_base = row['ref_base'].upper()
    
    # Skip lines below coverage threshold.
    if normal_depth < min_depth or tumour_depth < min_depth:
        return None
            
    normal_bases = get_bases(
                              ref_base,
                              row['normal_call_string'],
                              row['normal_base_qual_string'],
                              min_qual
                              )
    
    tumour_bases = get_bases(
                              ref_base,
                              row['tumour_call_string'],
                              row['tumour_base_qual_string'],
                              min_qual 
                              )
    
    tumour_non_ref_base, tumour_counts = get_counts(ref_base, tumour_bases)        
    normal_non_ref_base, normal_counts = get_counts(ref_base, normal_bases, non_ref_base=tumour_non_ref_base)        
    
    # Check again for lines below read depth. The first check above speeds things up, though redundant.
    d_N = normal_counts[0] + normal_counts[1]
    d_T = tumour_counts[0] + tumour_counts[1]

    if d_N < min_depth or d_T < min_depth:
        return None
    
    jcnt_entry = [ chr_coord, ref_base, normal_non_ref_base, tumour_non_ref_base ]
    jcnt_entry.extend(normal_counts)
    jcnt_entry.extend(tumour_counts)
    
    return jcnt_entry

def get_bases(ref_base, call_string, qual_string, min_qual):
    bases = parse_call_string(ref_base, call_string)
    
    quals = np.fromstring(qual_string, dtype=np.byte) - ascii_offset
    
    bases = np.array(bases)
    
    bases = bases[quals >= min_qual]
    
    return bases

def get_counts(ref_base, bases, non_ref_base=None):
    counter = Counter(bases)
    
    non_ref_base, counts = parse_counts(ref_base, counter, non_ref_base)
    
    return non_ref_base, counts

def parse_counts(ref_base, counter, non_ref_base=None):
    ref_counts = counter[ref_base]
    
    del counter[ref_base]
    del counter['N']
    
    # Check if there is any non-ref bases.
    if non_ref_base is not None:
        non_ref_counts = counter[non_ref_base]
    else:
        if len(counter) > 0:
            non_ref_base, non_ref_counts = counter.most_common(1)[0]
        else:
            non_ref_base = 'N'
            non_ref_counts = 0
    
    counts = (ref_counts, non_ref_counts)
    
    return non_ref_base, counts

def write_rows(jcnt_file, rows):
    for chr_name, chr_rows in rows.items():
        jcnt_file.add_rows(chr_name, chr_rows)
