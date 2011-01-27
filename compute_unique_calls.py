import csv
import glob
import os
from joint_snv_mix.file_formats.jsm import JointSnvMixReader

import urllib

def main( dir ):
    out_file_name = os.path.join( dir, 'out.tsv' )
    
    glob_str = os.path.join( dir, '*.tsv' )

    file_list = glob.glob( glob_str )
    
    sites = {}
    
    for file_name in file_list:
        if file_name == out_file_name:
            continue
        
        base_name = os.path.basename( file_name )
        case, tech, model, extension = base_name.split( '.' )
        
        print case, tech, model
        
        file_sites = load_sites_from_file( file_name, case )
        
        jsm_file_name = case + "." + tech + ".jsm"        
        jsm_file_name = os.path.join( dir, jsm_file_name )        
        jsm_sites = load_site_from_jsm( file_sites, jsm_file_name )
        
        for site in file_sites:
            if site not in sites:
                sites[site] = {}
                
                sites[site]['ref_base'] = jsm_sites[site]['ref_base']
                sites[site]['var_base'] = jsm_sites[site]['var_base']
                                
                sites[site]['method'] = set()
                
                sites[site]['excap_varscan_score'] = 'NA'
                sites[site]['excap_indep_score'] = 'NA'
                sites[site]['excap_joint_score'] = 'NA'
                
                sites[site]['solid_varscan_score'] = 'NA'
                sites[site]['solid_indep_score'] = 'NA'
                sites[site]['solid_joint_score'] = 'NA'
                
                sites[site]['excap_counts'] = ( 0, 0, 0, 0 )
                sites[site]['solid_counts'] = ( 0, 0, 0, 0 )
                
                chrom = site[1]
                position = site[2]                
                sites[site]['gene'] = get_gene_id( chrom, position, sites[site]['ref_base'], sites[site]['var_base'] )
            
            method = "_".join( ( tech, model ) )
            sites[site]['method'].add( method )
            
            score_type = "_".join( ( method, "score" ) )
            sites[site][score_type] = file_sites[site]
            
            count_str = "_".join( ( tech, 'counts' ) )
            sites[site][count_str] = jsm_sites[site]['counts']
            
            print sites[site]['gene'], case, model, tech
    
    rows = []
    for site in sorted( sites ):
        sites[site]['site'] = ",".join( [str( x ) for x in site] )
        sites[site]['name'] = '_'.join( [str( x ) for x in site] )
        
        sites[site]['excap_counts'] = ",".join( [str( x ) for x in sites[site]['excap_counts']] )
        sites[site]['solid_counts'] = ",".join( [str( x ) for x in sites[site]['solid_counts']] )
        
        sites[site]['method'] = ",".join( sorted( sites[site]['method'] ) )
        
        if len( sites[site]['gene'] ) is None:
            print sites[site]
            continue
        
        rows.append( sites[site] )

    fields = [
              'name',
              'site',
              'ref_base',
              'var_base',
              'gene',
              'excap_counts',
              'solid_counts',
              'method',
              'excap_varscan_score',
              'excap_indep_score',
              'excap_joint_score',
              'solid_varscan_score',
              'solid_indep_score',
              'solid_joint_score'
              ]
    
#    fields = [
#          'name',
#          'site',
#          'score',
#          'tech',
#          'model'
#          ]
    
    writer = csv.DictWriter( open( out_file_name, 'w' ), fields, delimiter='\t' )
    
    writer.writeheader()
    
    writer.writerows( rows )

def load_sites_from_file( file_name, case ):
    reader = csv.reader( open( file_name ), delimiter='\t' )
    
    sites = {}        

    for row in reader:
        chrom = row[0]
        position = int( row[1] )
        
#        site = "_".join( ( case, chrom, position ) )
        site = ( case, chrom, position )
        
        score = row[2]
        
        sites[site] = score
        
    return sites

def get_gene_id( chrom, pos, ref_base, var_base ):
    pos = str( pos )
    
    input = ",".join( ( chrom, pos, ref_base, var_base ) )

    url = "http://mutationassessor.org/?cm=var&var=" + input + "&fts=all&frm=txt"
    
    data = urllib.urlopen( url ).read()
    
    lines = data.split( '\n' )

    fields = lines[0].split( '\t' )
    values = lines[1].split( '\t' )
    
    processed_data = dict( zip( fields, values ) )
    
    if 'Non-coding' in fields:
        return None
    
    if 'Synonymous' in fields:
        return None
    
    if processed_data['Gene'] is None:
        raise Exception()
    
    return processed_data['Gene']

def load_site_from_jsm( search_sites, jsm_file_name ):
    jsm_reader = JointSnvMixReader( jsm_file_name )
    
    case = search_sites.keys()[0][0]
    
    site_set = set( [( x[1], x[2] ) for x in search_sites.keys()] )
    
    chroms = jsm_reader.get_chr_list()
    
    found_sites = {}
    
    for chrom in chroms:        
        rows = jsm_reader.get_rows( chrom )
        
        for row in rows:        
            pos = row[0]
            
            if ( chrom, pos ) in site_set:
                site = ( case, chrom, pos )
                
                found_sites[site] = {}
                found_sites[site]['ref_base'] = row[1]
                found_sites[site]['var_base'] = row[3]
                found_sites[site]['counts'] = ( row[4], row[5], row[6], row[7] )
    
    print set( search_sites.keys() ) - set( found_sites.keys() )

    jsm_reader.close()
                
    return found_sites

        

if __name__ == "__main__":
    import sys
    
    dir = sys.argv[1]
    
    main( dir )
