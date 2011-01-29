import csv
import glob
import os
from joint_snv_mix.file_formats.jsm import JointSnvMixReader

import urllib
import multiprocessing

def main( dir ):
    out_file_name = os.path.join( dir, 'out_p.tsv' )
    
    glob_str = os.path.join( dir, '*.tsv' )

    file_list = glob.glob( glob_str )
    
    if 'out.tsv' in file_list:
        file_list.remove( 'out.tsv' )
    
    if out_file_name in file_list:
        file_list.remove( out_file_name )
    
    p = multiprocessing.Pool()    
    results = p.map( get_sites, file_list )
    p.close()
    
    all_sites = merge_sites( results )
    
    rows = []
    
    for site in sorted( all_sites ):
        all_sites[site]['site'] = ",".join( [str( x ) for x in site] )
        all_sites[site]['name'] = '_'.join( [str( x ) for x in site] )
        
        all_sites[site]['excap_counts'] = ",".join( [str( x ) for x in all_sites[site]['excap_counts']] )
        all_sites[site]['solid_counts'] = ",".join( [str( x ) for x in all_sites[site]['solid_counts']] )
        
        all_sites[site]['method'] = ",".join( sorted( all_sites[site]['method'] ) )
        
        if len( all_sites[site]['gene'] ) is None:
            print all_sites[site]
            continue
        
        rows.append( all_sites[site] )

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
              'solid_joint_score',
              'excap_joint_probs',
              'solid_joint_probs'              
              ]
    
    writer = csv.DictWriter( open( out_file_name, 'w' ), fields, delimiter='\t' )
    
    writer.writeheader()
    
    writer.writerows( rows )
    
def get_sites( file_name ):
    sites = {}
    
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
            
            sites[site]['excap_counts'] = "NA"
            sites[site]['solid_counts'] = "NA"
            
            sites[site]['excap_joint_probs'] = "NA"
            sites[site]['solid_joint_probs'] = "NA"
            
            chrom = site[1]
            position = site[2]                
            sites[site]['gene'] = get_gene_id( chrom, position, sites[site]['ref_base'], sites[site]['var_base'] )
        
        method = "_".join( ( tech, model ) )
        sites[site]['method'].add( method )
        
        score_type = "_".join( ( method, "score" ) )
        sites[site][score_type] = file_sites[site]
        
        count_str = "_".join( ( tech, 'counts' ) )
        sites[site][count_str] = jsm_sites[site]['counts']
        
        jsm_score_str = "_".join( ( tech, "joint_probs" ) )
        sites[site][jsm_score_str] = jsm_sites[site]['probs']
                    
        print sites[site]['gene'], case, model, tech
        
    return sites

def merge_sites( sites ):
    merged_sites = {}
    
    for case_sites in sites:
        for site in case_sites: 
                       
            if site not in merged_sites:
                merged_sites[site] = {}
            
            for row_name in case_sites[site]:
                merged_sites[site][row_name] = case_sites[site][row_name]

    return merged_sites

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
                found_sites[site]['probs'] = ( row[8], row[9], row[10],
                                               row[11], row[12], row[13],
                                               row[14], row[15], row[16] )
    
    print set( search_sites.keys() ) - set( found_sites.keys() )

    jsm_reader.close()
                
    return found_sites

        

if __name__ == "__main__":
    import sys
    
    dir = sys.argv[1]
    
    main( dir )
