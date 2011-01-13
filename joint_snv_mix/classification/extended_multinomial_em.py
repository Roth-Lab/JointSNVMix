'''
Created on 2010-12-09

@author: Andrew Roth
'''
import ConfigParser

import numpy as np

from joint_snv_mix.file_formats.mcnt import MultinomialCountsReader
from joint_snv_mix import constants

from jsm_models.em.extended_multinomial.extended_multinomial import JointExtendedMultinomialModel
from jsm_models.utils.subsample import subsample
from jsm_models.data import MultinomialData
from joint_snv_mix.file_formats.jemm import JointExtendedMultiMixWriter

nclass = 15

def main( args ):
    model = JointExtendedMultinomialModel()
    
    reader = MultinomialCountsReader( args.jcnt_file_name )
    writer = JointExtendedMultiMixWriter( args.jsm_file_name )
    
    # Load parameters by training or from file.
    if args.train:
        if args.subsample_size > 0:
            counts = subsample( reader, args.subsample_size )
        else:
            counts = reader.get_counts()
        
        data = MultinomialData( counts )
        
        priors = parse_priors_file( args.priors_file_name, data.nrows )
        parameters = model.train( data, priors, args.max_iters, args.convergence_threshold )   
    else:
        parameters = parse_parameters_file( args.params_file_name )
    
    chr_list = reader.get_chr_list()
    
    for chr_name in sorted( chr_list ):
        counts = reader.get_counts( chr_name )
        
        data = MultinomialData( counts )
        
        responsibilities = model.classify( data, parameters )
        
        jcnt_rows = reader.get_rows( chr_name )
        
        writer.write_data( chr_name, jcnt_rows, responsibilities )
        
    reader.close()
    writer.close()
        
def parse_priors_file( priors_file_name, n ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( priors_file_name )
    
    priors = {}
    priors['kappa'] = np.zeros( ( nclass ** 2, ) )
    priors['delta'] = np.zeros( ( 2, nclass, 4 ) )
    
    for i, genotype_tuple in enumerate( constants.joint_extended_multinomial_genotypes ):
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
        
        
    for i, genotype in enumerate( constants.extended_multinomial_genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        for j, nuc in enumerate( constants.nucleotides ):
            normal_genotype_nuc = "_".join( ( normal_genotype, nuc ) )
            tumour_genotype_nuc = "_".join( ( tumour_genotype, nuc ) )
        
            priors['delta'][0, i, j] = parser.getfloat( 'delta', normal_genotype_nuc )
            priors['delta'][1, i, j] = parser.getfloat( 'delta', tumour_genotype_nuc )
    
    priors['delta'] = np.log( n ) * priors['delta']
    
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
