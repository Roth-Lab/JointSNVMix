#!/usr/bin/env python

import csv

import pysam

from bam_counter.joint_binary_counter import JointBinaryBaseCounter
from bam_counter.positions_iterator import PositionCountIterator, CharacterDelimitedPositionsReader

from joint_snv_mix.file_formats.jcnt import JointCountsWriter

def paired_bams_to_jcnt(args):
    normal_bam = pysam.Samfile(args.normal_bam_file_name, 'rb')
    tumour_bam = pysam.Samfile(args.tumour_bam_file_name, 'rb')

    ref_genome_fasta = pysam.Fastafile(args.reference_genome_file_name)
    
    counter = JointBinaryBaseCounter(normal_bam, tumour_bam, ref_genome_fasta, args.min_base_qual, args.min_map_qual)
    
    writer = JointCountsWriter(args.jcnt_file_name)
    
    if args.positions_file is None:
        write_all_positions(counter, writer, args.min_depth)
    else:
        write_specified_positions(counter, writer, args.min_depth, args.positions_file, args.delimiter)
    
    writer.close()
        
def write_all_positions(counter, writer, min_depth):
    for ref in sorted(counter.refs):
        print ref
        for row in counter.iter_ref(ref):
            if row.depth >= min_depth:        
                out_row = format_row(row)                
                writer.add_row(row.ref, out_row)

def write_specified_positions(counter, writer, min_depth, positions_file, delimiter):
    if delimiter == 'tab':
        delimiter_char = '\t'
    elif delimiter == 'space':
        delimiter_char = ' '
    elif delimiter == 'comma':
        delimiter_char = ','
        
    positions_reader = CharacterDelimitedPositionsReader(positions_file, delimiter_char)
    
    positions_iter = PositionsIterator(counter, positions_reader)
    
    for row in positions_iter:
        if row.depth >= min_depth:
            writer.write(str(row))
            writer.write("\n")
            
def format_row(row):
    out_row = [row.position, row.ref_base, row.non_ref_base]
    out_row.extend(row.counts)    
    
    return out_row