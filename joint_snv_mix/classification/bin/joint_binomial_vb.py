'''
Created on 2010-12-09

@author: Andrew Roth
'''
import ConfigParser

import numpy as np

from joint_snv_mix import constants
from joint_snv_mix.classification.data import JointData
from joint_snv_mix.classification.vb.joint_models.binomial.mixture_model import MixtureModel as JointBinomialModel
from joint_snv_mix.classification.utils.subsample import subsample
from joint_snv_mix.file_formats.jcnt import JointCountsReader
from joint_snv_mix.file_formats.jsm import JointSnvMixWriter


def main( args ):
    model = JointBinomialModel()
    
    reader = JointCountsReader( args.jcnt_file_name )
    writer = JointSnvMixWriter( args.jsm_file_name )
    
    # Load parameters by training or from file.
    if args.train:
        if args.subsample_size > 0:
            counts = subsample( reader, args.subsample_size )
        else:
            counts = reader.get_counts()
        
        data = JointData( counts )
        
        priors = parse_priors_file( args.priors_file, data.nrows )
        posterior = model.train( data, priors, args.max_iters, args.convergence_threshold )   
    else:
        parameters = parse_parameters_file( args.params_file_name )
    
    chr_list = reader.get_chr_list()
    
    for chr_name in sorted( chr_list ):
        counts = reader.get_counts( chr_name )
        
        data = JointData( counts )
        
        responsibilities = model.classify( data, posterior )
        
        jcnt_rows = reader.get_rows( chr_name )
        
        writer.write_data( chr_name, jcnt_rows, responsibilities )
        
    reader.close()
    writer.close()
        
def parse_priors_file( priors_file_name, n ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( priors_file_name )
    
    priors = {}
    priors['kappa'] = np.zeros( ( 9, ) )
    priors['alpha'] = np.zeros( ( 2, 3 ) )
    priors['beta'] = np.zeros( ( 2, 3 ) )
    
    for i, genotype_tuple in enumerate( constants.joint_genotypes ):
        genotype = "_".join( genotype_tuple )
        
        priors['kappa'][i] = parser.getfloat( 'kappa', genotype )
        
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
        
        priors['alpha'][0, i] = parser.getfloat( 'alpha', normal_genotype )
        priors['beta'][0, i] = parser.getfloat( 'beta', normal_genotype )
        
        priors['alpha'][1, i] = parser.getfloat( 'alpha', tumour_genotype )
        priors['beta'][1, i] = parser.getfloat( 'beta', tumour_genotype )
        
    return priors
    
def parse_parameters_file( parameters_file_name ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( parameters_file_name )
    
    parameters = {}
    parameters['mu'] = np.zeros( ( 2, 3 ) )
    parameters['pi'] = np.zeros( ( 9, ) )
    
    for i, genotype_tuple in enumerate( constants.joint_genotypes ):
        genotype = "_".join( genotype_tuple )
        
        parameters['pi'][i] = parser.getfloat( 'pi', genotype )
    
    parameters['pi'] = parameters['pi'] / parameters['pi'].sum()
        
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        parameters['mu'][0, i] = parser.getfloat( 'mu', normal_genotype )
        parameters['mu'][1, i] = parser.getfloat( 'mu', tumour_genotype )
        
    return parameters
