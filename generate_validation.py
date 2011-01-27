import csv
import glob
import os
from joint_snv_mix.file_formats.jsm import JointSnvMixReader
from joint_snv_mix.classification.classification import run_classifier
from argparse import Namespace

excluded_chrom = ['MT', 'Y']

def main( validation_dir ):    
    varscan_predictions = load_varscan_predictions( validation_dir )
    
    jsm_files, ism_files = run_jsm( validation_dir )
    
    joint_predictions = load_jsm_predictions( jsm_files, varscan_predictions )
    
    indep_predictions = load_jsm_predictions( ism_files, varscan_predictions )
    
    compare_positions( varscan_predictions, validation_dir, "varscan" )
    compare_positions( indep_predictions, validation_dir, "indep" )
    compare_positions( joint_predictions, validation_dir, "joint" )
    
def load_varscan_predictions( validation_dir ):
    glob_str = os.path.join( validation_dir, "*.hc" )
    
    varscan_file_list = glob.glob( glob_str )
    
    predictions = {}
    
    for file_name in sorted( varscan_file_list ):        
        base_name = os.path.basename( file_name )
        
        case = base_name.split( '.' )[0]
        
        if case not in predictions.keys():
            predictions[case] = {}
            
        tech = base_name.split( '.' )[1]
        
        predictions[case][tech] = load_predictions_from_varscan_file( file_name )
        
    return predictions
        
def load_predictions_from_varscan_file( file_name ):
    reader = csv.DictReader( open( file_name ), delimiter='\t' )
    
    positions = []
    
    for row in reader:
        if row['chrom'] in excluded_chrom:
            continue
        
        positions.append( ( 
                          row['chrom'],
                          int( row['position'] ),
                          float( row['somatic_p_value'] )
                          ) )
        
    positions = sorted( positions, key=lambda position: position[2] )
    
    return positions

def run_jsm( validation_dir ):
    glob_str = os.path.join( validation_dir, "*.jcnt" )
    
    jcnt_file_list = glob.glob( glob_str )
    
    jsm_files = {}
    
    ism_files = {}
    
    for jcnt_file_name in jcnt_file_list:
        base_name = os.path.basename( jcnt_file_name )
        
        case = base_name.split( '.' )[0]
        
        if case not in jsm_files.keys():
            jsm_files[case] = {}
            ism_files[case] = {}
            
        tech = base_name.split( '.' )[1]
        
        jsm_file_name = case + "." + tech + ".jsm"
        jsm_file_name = os.path.join( validation_dir, jsm_file_name )
               
        jsm_files[case][tech] = jsm_file_name
        
        if not os.path.exists( jsm_file_name ):
            run_joint_bb( jcnt_file_name, jsm_file_name )
                
        ism_file_name = case + "." + tech + ".ism"
        ism_file_name = os.path.join( validation_dir, ism_file_name )
        
        ism_files[case][tech] = ism_file_name
        
        if not os.path.exists( ism_file_name ):
            run_indep_bin( jcnt_file_name, ism_file_name )        
    
    return jsm_files, ism_files

def run_joint_bb( jcnt_file_name, jsm_file_name ):
    args = {}
    
    args['model'] = "joint"
    args['density'] = "beta_binomial"
    args['priors_file'] = '/home/andrew/workspace/joint_snv_mix/config/joint_bb_em.priors.cfg'
    
    args['jcnt_file_name'] = jcnt_file_name
    args['jsm_file_name'] = jsm_file_name
    
    args['max_iters'] = int( 1e3 ) 
    args['convergence_threshold'] = float( 1e-6 )
    
    args['train'] = True
    args['subsample_size'] = int( 1e6 )
    
    args = Namespace( **args )
    
    run_classifier( args )

def run_indep_bin( jcnt_file_name, ism_file_name ):
    args = {}
    
    args['model'] = "independent"
    args['density'] = "binomial"
    args['priors_file'] = '/home/andrew/workspace/joint_snv_mix/config/indep_bin_em.priors.cfg'
    
    args['jcnt_file_name'] = jcnt_file_name
    args['jsm_file_name'] = ism_file_name
    
    args['max_iters'] = int( 1e3 ) 
    args['convergence_threshold'] = float( 1e-6 )
    
    args['train'] = True
    args['subsample_size'] = int( 1e6 )
    
    args = Namespace( **args )
    
    run_classifier( args )
    
def load_jsm_predictions( jsm_files, varscan_predictions ):
    jsm_predictions = {}
    
    for case in jsm_files.keys():
        jsm_predictions[case] = {}
        
        for tech in jsm_files[case].keys():
            print "Loading ", case, tech
            
            n = len( varscan_predictions[case][tech] )
            
            predictions = load_predictions_from_jsm_file( jsm_files[case][tech] )
            
            predictions = predictions[:n]
            
            jsm_predictions[case][tech] = predictions
    
    return jsm_predictions
    

def load_predictions_from_jsm_file( jsm_file_name ):
    position_score = load_somatics( jsm_file_name )
    
    position_score = sort_position_score( position_score )    
    
    return position_score

def load_somatics( jsm_file_name ):
    reader = JointSnvMixReader( jsm_file_name )
    
    chr_list = reader.get_chr_list()
    
    position_score = []
    
    for chr_name in sorted( chr_list ):
        if chr_name in excluded_chrom:
            continue
#        
        chr_rows = reader.get_rows( chr_name )
        
        for row in chr_rows:
            position = int( row['position'] )
            score = row['p_aa_ab'] + row['p_aa_bb']
            
            position_score.append( ( chr_name, position, score ) )
    
    reader.close()
        
    return position_score

def sort_position_score( position_score ):
    position_score = sorted( position_score, key=lambda x: x[2], reverse=True )
    
    return position_score

def compare_positions( predictions, dir, method_name ):
    for case in predictions.keys():
        for tech in predictions[case].keys():
            file_name = case + "." + tech + "." + method_name + ".tsv"
            
            file_name = os.path.join( dir, file_name )
        
            writer = csv.writer( open( file_name, 'w' ), delimiter='\t' )
            
            writer.writerows( predictions[case][tech] )


def get_position_set( predictions ):
    positions = set()
    
    for prediction in predictions:
        positions.add( tuple( prediction[:2] ) )
    
    return positions
        
    
if __name__ == "__main__":
    import sys
    
    validation_dir = sys.argv[1]
    
    main( validation_dir )
