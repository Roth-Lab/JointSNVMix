import csv
import glob
import os
from joint_snv_mix.file_formats.jsm import JointSnvMixReader
from joint_snv_mix.classification.classification import run_classifier
from argparse import Namespace
import bisect
import multiprocessing

excluded_chrom = ['MT', 'Y']

def main( validation_dir ):    
    varscan_predictions = load_varscan_predictions( validation_dir )
    write_predictions_to_tsv( varscan_predictions, validation_dir, "varscan" )
    
    models = {
              'indep_bin' : run_indep_bin,
              'joint_bin' : run_joint_bin,
              'joint_bb' : run_joint_bb,
              'joint_bb_fixed' : run_joint_bb_fixed              
              }
    
    p = multiprocessing.Pool()
    
    for model_name, run_func in models.iteritems():
        args = ( model_name, run_func, varscan_predictions )
        p.apply_async( write_model_predicitions, args )

def write_model_predicitions( model_name, run_func, varscan_predictions ):
    files = run_model( validation_dir, model_name, run_func )
    predictions = load_jsm_predictions( files, varscan_predictions )
    write_predictions_to_tsv( predictions, validation_dir, model_name )
            
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

def run_model( validation_dir, model_name, run_model_func ):
    glob_str = os.path.join( validation_dir, "*.jcnt" )
    
    jcnt_file_list = glob.glob( glob_str )
    
    model_files = {}
    
    for jcnt_file_name in jcnt_file_list:
        base_name = os.path.basename( jcnt_file_name )
        
        case = base_name.split( '.' )[0]
        tech = base_name.split( '.' )[1]
        
        if case not in model_files.keys():
            model_files[case] = {}    

        model_file_name = case + "." + tech + "." + model_name + ".jsm" 
        
        model_file_name = os.path.join( validation_dir, model_file_name )
        
        model_files[case][tech] = model_file_name
        
        if not os.path.exists( model_file_name ):
            try:
                run_model_func( jcnt_file_name, model_file_name )
            except:
                pass 

def run_indep_bin( jcnt_file_name, jsm_file_name ):
    args = {}
    
    args['model'] = "independent"
    args['density'] = "binomial"
    args['priors_file'] = '/home/andrew/workspace/joint_snv_mix/config/indep_bin.priors.cfg'
    
    args['jcnt_file_name'] = jcnt_file_name
    args['jsm_file_name'] = jsm_file_name
    
    args['max_iters'] = int( 1e3 ) 
    args['convergence_threshold'] = float( 1e-6 )
    
    args['train'] = True
    args['subsample_size'] = int( 1e6 )
    
    args = Namespace( **args )
    
    run_classifier( args )
    
def run_indep_bb( jcnt_file_name, jsm_file_name ):
    args = {}
    
    args['model'] = "independent"
    args['density'] = "beta_binomial"
    args['priors_file'] = '/home/andrew/workspace/joint_snv_mix/config/indep_bb.priors.cfg'
    
    args['jcnt_file_name'] = jcnt_file_name
    args['jsm_file_name'] = jsm_file_name
    
    args['max_iters'] = int( 1e3 ) 
    args['convergence_threshold'] = float( 1e-6 )
    
    args['train'] = True
    args['subsample_size'] = int( 1e6 )
    
    args = Namespace( **args )
    
    run_classifier( args )
    
def run_joint_bin( jcnt_file_name, jsm_file_name ):
    args = {}
    
    args['model'] = "joint"
    args['density'] = "binomial"
    args['priors_file'] = '/home/andrew/workspace/joint_snv_mix/config/joint_bin.priors.cfg'
    
    args['jcnt_file_name'] = jcnt_file_name
    args['jsm_file_name'] = jsm_file_name
    
    args['max_iters'] = int( 1e3 ) 
    args['convergence_threshold'] = float( 1e-6 )
    
    args['train'] = True
    args['subsample_size'] = int( 1e6 )
    
    args = Namespace( **args )
    
    run_classifier( args )

def run_joint_bb( jcnt_file_name, jsm_file_name ):
    args = {}
    
    args['model'] = "joint"
    args['density'] = "beta_binomial"
    args['priors_file'] = '/home/andrew/workspace/joint_snv_mix/config/joint_bb.priors.cfg'
    
    args['jcnt_file_name'] = jcnt_file_name
    args['jsm_file_name'] = jsm_file_name
    
    args['max_iters'] = int( 1e3 ) 
    args['convergence_threshold'] = float( 1e-6 )
    
    args['train'] = True
    args['subsample_size'] = int( 1e6 )
    
    args = Namespace( **args )
    
    run_classifier( args )

def run_joint_bb_fixed( jcnt_file_name, jsm_file_name ):
    args = {}
    
    args['model'] = "joint"
    args['density'] = "beta_binomial"
    args['params_file'] = '/home/andrew/workspace/joint_snv_mix/config/joint_bb.params.cfg'
    
    args['jcnt_file_name'] = jcnt_file_name
    args['jsm_file_name'] = jsm_file_name
    
    args['train'] = False
    
    args = Namespace( **args )
    
    run_classifier( args )
    
def load_jsm_predictions( jsm_files, varscan_predictions ):
    jsm_predictions = {}
    
    for case in jsm_files.keys():
        jsm_predictions[case] = {}
        
        for tech in jsm_files[case].keys():
            print "Loading ", case, tech
            
            n = len( varscan_predictions[case][tech] )
            
            predictions = load_somatics( jsm_files[case][tech], n )
            
            jsm_predictions[case][tech] = predictions
    
    return jsm_predictions

def load_somatics( jsm_file_name, n ):
    reader = JointSnvMixReader( jsm_file_name )
    
    chr_list = reader.get_chr_list()
    
    position_score = []
    scores = []
    
    for chr_name in sorted( chr_list ):
        if chr_name in excluded_chrom:
            continue
                
        chr_rows = reader.get_rows( chr_name )
        
        for row in chr_rows:
            position = int( row['position'] )
            score = float( row['p_aa_ab'] + row['p_aa_bb'] )            
            
            insert_position = bisect.bisect( scores, score )
            
            if insert_position > 0 or len( scores ) == 0:               
                position_score.insert( insert_position, ( chr_name, position, score ) )
                scores.insert( insert_position, score )
                
                if len( scores ) >= n:
                    scores.pop( 0 )
                    position_score.pop( 0 )                        
                        
    reader.close()
        
    return position_score

def write_predictions_to_tsv( predictions, dir, method_name ):
    for case in predictions.keys():
        for tech in predictions[case].keys():
            file_name = case + "." + tech + "." + method_name + ".tsv"
            
            file_name = os.path.join( dir, file_name )
        
            writer = csv.writer( open( file_name, 'w' ), delimiter='\t' )
            
            writer.writerows( predictions[case][tech] )
        
    
if __name__ == "__main__":
    import sys
    
    validation_dir = sys.argv[1]
    
    main( validation_dir )
