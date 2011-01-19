'''
Created on 2010-12-09

@author: Andrew Roth
'''
import ConfigParser

import numpy as np

from joint_snv_mix import constants
from joint_snv_mix.classification.data import IndependentData
from joint_snv_mix.classification.em.independent_models.beta_binomial import IndependenBetaBinomialModel
from joint_snv_mix.classification.utils.subsample import subsample
from joint_snv_mix.file_formats.jcnt import JointCountsReader
from joint_snv_mix.file_formats.jsm import JointSnvMixWriter

def main( args ):
    model = IndependenBetaBinomialModel()
    
    reader = JointCountsReader( args.jcnt_file_name )
    writer = JointSnvMixWriter( args.jsm_file_name )
    
    # Load parameters by training or from file.
    if args.train:
        if args.subsample_size > 0:
            counts = subsample( reader, args.subsample_size )
        else:
            counts = reader.get_counts()
        
        normal_data = IndependentData( counts, 'normal' )
        tumour_data = IndependentData( counts, 'tumour' )
        
        n = normal_data.nrows 
        
        normal_priors, tumour_priors = parse_priors_file( args.priors_file, n )
        
        normal_parameters = model.train( normal_data, normal_priors, args.max_iters, args.convergence_threshold )
        tumour_parameters = model.train( tumour_data, tumour_priors, args.max_iters, args.convergence_threshold )   
    else:
        normal_parameters, tumour_data = parse_parameters_file( args.params_file_name )
    
    chr_list = reader.get_chr_list()
    
    for chr_name in sorted( chr_list ):
        counts = reader.get_counts( chr_name )
        
        normal_data = IndependentData( counts, 'normal' )
        tumour_data = IndependentData( counts, 'tumour' )
        
        normal_responsibilities = model.classify( normal_data, normal_parameters )
        tumour_responsibilities = model.classify( tumour_data, tumour_parameters )
        
        responsibilities = get_joint_responsibilities( normal_responsibilities, tumour_responsibilities )
              
        jcnt_rows = reader.get_rows( chr_name )
        
        writer.write_data( chr_name, jcnt_rows, responsibilities )
        
    reader.close()
    writer.close()
        
def parse_priors_file( priors_file_name, n ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( priors_file_name )
    
    priors = {}
    priors['kappa'] = np.zeros( ( 2, 3 ) )
    priors['location'] = np.zeros( ( 2, 3, 2 ) )
    priors['precision'] = np.zeros( ( 2, 3, 2 ) )
    
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        priors['kappa'][0, i] = parser.getfloat( 'kappa', normal_genotype )
        priors['kappa'][1, i] = parser.getfloat( 'kappa', tumour_genotype )
        
    scaling = parser.getint( 'kappa', 'scaling' )
    
    if scaling == 1:
        priors['kappa'] = np.log( n ) * priors['kappa']
    elif scaling == 2:
        priors['kappa'] = np.sqrt( n ) * priors['kappa']
    elif scaling == 3:
        priors['kappa'] = n / 10. * priors['kappa']
    elif scaling == 4:
        priors['kappa'] = n / 2. * priors['kappa']
    elif scaling == 4:
        priors['kappa'] = n * priors['kappa']
        
        
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        priors['location'][0, i, 0] = parser.getfloat( 'location_alpha', normal_genotype )
        priors['location'][0, i, 1] = parser.getfloat( 'location_beta', normal_genotype )
        
        priors['location'][1, i, 0] = parser.getfloat( 'location_alpha', tumour_genotype )
        priors['location'][1, i, 1] = parser.getfloat( 'location_beta', tumour_genotype )
        
        priors['precision'][0, i, 0] = parser.getfloat( 'precision_shape', normal_genotype )
        priors['precision'][0, i, 1] = parser.getfloat( 'precision_scale', normal_genotype )
        
        priors['precision'][1, i, 0] = parser.getfloat( 'precision_shape', tumour_genotype )
        priors['precision'][1, i, 1] = parser.getfloat( 'precision_scale', tumour_genotype )
        
    normal_priors = {'kappa' : priors['kappa'][0], 'location' : priors['location'][0], 'precision' : priors['precision'][0]}
    tumour_priors = {'kappa' : priors['kappa'][1], 'location' : priors['location'][1], 'precision' : priors['precision'][1]}
        
    return normal_priors, tumour_priors
    
def parse_parameters_file( parameters_file_name ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( parameters_file_name )
    
    parameters = {}
    parameters['alpha'] = np.zeros( ( 2, 3 ) )
    parameters['beta'] = np.zeros( ( 2, 3 ) )
    parameters['pi'] = np.zeros( ( 2, 3 ) )
    
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        parameters['pi'][0, i] = parser.getfloat( 'pi', normal_genotype )
        parameters['pi'][1, i] = parser.getfloat( 'pi', tumour_genotype )
    
    parameters['pi'][0, :] = parameters['pi'][0, :] / parameters['pi'][0, :].sum()
    parameters['pi'][1, :] = parameters['pi'][1, :] / parameters['pi'][1, :].sum()
        
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        parameters['alpha'][0, i] = parser.getfloat( 'alpha', normal_genotype )
        parameters['alpha'][1, i] = parser.getfloat( 'alpha', tumour_genotype )
        
        parameters['beta'][0, i] = parser.getfloat( 'beta', normal_genotype )
        parameters['beta'][1, i] = parser.getfloat( 'beta', tumour_genotype )
        
    normal_parameters = {'pi' : parameters['pi'][0], 'alpha' : parameters['alpha'][0], 'beta' : parameters['beta'][0] }
    tumour_parameters = {'pi' : parameters['pi'][1], 'alpha' : parameters['alpha'][1], 'beta' : parameters['beta'][1] }
        
    return normal_parameters, tumour_parameters

def get_joint_responsibilities( normal_responsibilities, tumour_responsibilities ):
    eps = np.finfo( np.double ).eps
    
    normal_responsibilities = np.log( normal_responsibilities )
    tumour_responsibilities = np.log( tumour_responsibilities )
      
    column_shape = ( normal_responsibilities[:, 0].size, 1 )
    
    log_responsibilities = np.hstack( ( 
                             normal_responsibilities[:, 0].reshape( column_shape ) + tumour_responsibilities ,
                             normal_responsibilities[:, 1].reshape( column_shape ) + tumour_responsibilities ,
                             normal_responsibilities[:, 2].reshape( column_shape ) + tumour_responsibilities
                             ) )
    
    responsibilities = np.exp( log_responsibilities )
    
    responsibilities[responsibilities <= eps] = 0.
    
    return responsibilities

#=======================================================================================================================
# Classes
#=======================================================================================================================

        


