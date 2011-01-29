#!/usr/bin/env python
import csv

from joint_snv_mix.file_formats.jsm import JointSnvMixReader

def main( args ):    
    positions = get_positions( args.positions_file_name )
    
    get_position_probabilities( args.jsm_file_name, positions )

def get_positions( position_file_name ):
    fileds = ['chr_name', 'coordinates']
    reader = csv.DictReader( open( position_file_name ), delimiter='\t', fieldnames=fileds )
    
    positions = {}
    
    for row in reader:
        chr_name = row['chr_name']
        coordinate = row['coordinates']        
        
        if chr_name not in positions:
            positions[chr_name] = []
        
        positions[chr_name].append( coordinate )
    
    return positions
    
def get_position_probabilities( jsm_file_name, positions ):
    reader = JointSnvMixReader( jsm_file_name )
    
    for chromosome, coordinates in sorted( positions.items() ):
        for coord in coordinates:
            try:
                jsm_data = reader.get_position( chromosome, coord )
            except KeyError:
                print chromosome, " not in jsm file."
                continue
            
            if jsm_data:
                row = list( jsm_data )
    
                row[0] = chromosome + ":" + str( row[0] )
                
                print "\t".join( [str( x ) for x in row] )
            else:
                print chromosome + ":" + coord, " not if jsm file."
        
    reader.close()
