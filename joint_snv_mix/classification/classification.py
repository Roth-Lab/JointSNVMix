'''
Created on 2011-01-13

@author: Andrew Roth
'''
import ConfigParser
import math
import random

import numpy as np

from joint_snv_mix import constants
from joint_snv_mix.file_formats.jcnt import JointCountsReader
from joint_snv_mix.file_formats.jsm import JointSnvMixWriter
from joint_snv_mix.classification.data import IndependentData, JointData
from joint_snv_mix.classification.models import IndependentBinomialModel, JointBinomialModel, JointBetaBinomialModel, \
    IndependenBetaBinomialModel

def run_classifier( args ):
    if args.model == "independent":
        run_independent_model( args )
            
    elif args.model == "joint":
        run_joint_model( args )

    elif args.model == "chromosome":
        run_chromosome_model( args )

def subsample( reader, sample_size ):
    chr_list = reader.get_chr_list()
    
    sample = []
    
    nrows = reader.get_data_set_size()
    
    for chr_name in chr_list:
        chr_size = reader.get_chr_size( chr_name=chr_name )
        
        chr_sample_size = math.ceil( float( chr_size ) / nrows * sample_size )
        
        chr_sample_size = int( chr_sample_size )
        
        chr_sample_size = min( chr_size, chr_sample_size )
        
        chr_sample_indices = random.sample( xrange( chr_size ), chr_sample_size )
        
        chr_counts = reader.get_counts( chr_name )
        
        chr_sample = chr_counts[chr_sample_indices]
        
        sample.append( chr_sample )
        
    sample = np.vstack( sample )
    
    return sample

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

def scale_priors( priors, n, scaling ):
    if scaling == 1:
        priors['kappa'] = np.log( n ) * priors['kappa']
    elif scaling == 2:
        priors['kappa'] = np.sqrt( n ) * priors['kappa']
    elif scaling == 3:
        priors['kappa'] = n / 10. * priors['kappa']
    elif scaling == 4:
        priors['kappa'] = n / 2. * priors['kappa']
    elif scaling == 5:
        priors['kappa'] = n * priors['kappa']
        
    return priors

#=======================================================================================================================
# Independent
#=======================================================================================================================

def run_independent_model( args ):
    if args.density == "beta_binomial":
        model = IndependenBetaBinomialModel()
    elif args.density == "binomial":
        model = IndependentBinomialModel()
    else:
        raise NotImplemented
    
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
        
        priors = parse_independent_priors_file( args.priors_file, n, args.density )
        
        normal_parameters = model.train( normal_data, priors['normal'], args.max_iters, args.convergence_threshold )
        tumour_parameters = model.train( tumour_data, priors['tumour'], args.max_iters, args.convergence_threshold )   
    else:
        normal_parameters, tumour_data = parse_independent_parameters_file( args.params_file, args.density )
    
    chr_list = reader.get_chr_list()
    
    for chr_name in sorted( chr_list ):
        counts = reader.get_counts( chr_name )
        jcnt_rows = reader.get_rows( chr_name )
        
        end = reader.get_chr_size( chr_name )

        n = int( 1e5 )
        start = 0
        stop = min( n, end )
        

        while start < end:
            sub_counts = counts[start:stop]
            sub_rows = jcnt_rows[start:stop]
                          
            normal_data = IndependentData( sub_counts, 'normal' )
            tumour_data = IndependentData( sub_counts, 'tumour' )
            
            normal_responsibilities = model.classify( normal_data, normal_parameters )
            tumour_responsibilities = model.classify( tumour_data, tumour_parameters )
            
            responsibilities = get_joint_responsibilities( normal_responsibilities, tumour_responsibilities )
        
            writer.write_data( chr_name, sub_rows, responsibilities )
            
            start = stop
            stop = min( stop + n, end )
        
    reader.close()
    writer.close()
        
def parse_independent_priors_file( priors_file_name, n, density ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( priors_file_name )
    
    priors = {}
    
    priors['normal'] = {}
    priors['tumour'] = {}
    
    priors['normal']['kappa'] = np.zeros( ( 3, ) )
    priors['tumour']['kappa'] = np.zeros( ( 3, ) )
    
    for genome in ( 'normal', 'tumour' ):
        for i, genotype in enumerate( constants.genotypes ):
            genome_genotype = "_".join( ( genome, genotype ) )
            
            priors[genome]['kappa'][i] = parser.getfloat( 'kappa', genome_genotype )
        
    scaling = parser.getint( 'kappa', 'scaling' )
    
    priors['normal'] = scale_priors( priors['normal'], n, scaling )
    priors['tumour'] = scale_priors( priors['tumour'], n, scaling )
    
    if density == "beta_binomial":
        priors = get_independent_beta_binomial_density_priors( parser, priors )
    elif density == "binomial":
        priors = get_independent_binomial_density_priors( parser, priors )
    else:
        raise NotImplementedError
    
    return priors

def get_independent_beta_binomial_density_priors( parser, priors ):
    shape = ( 3, )
    
    priors['normal']['precision'] = {}
    priors['normal']['precision']['shape'] = np.zeros( shape )
    priors['normal']['precision']['scale'] = np.zeros( shape )
    
    priors['tumour']['precision'] = {}
    priors['tumour']['precision']['shape'] = np.zeros( shape )
    priors['tumour']['precision']['scale'] = np.zeros( shape )
    
    
    priors['normal']['location'] = {}
    priors['normal']['location']['alpha'] = np.zeros( shape )
    priors['normal']['location']['beta'] = np.zeros( shape )
    
    priors['tumour']['location'] = {}
    priors['tumour']['location']['alpha'] = np.zeros( shape )
    priors['tumour']['location']['beta'] = np.zeros( shape )
    
    for genome in ( 'normal', 'tumour' ):
        for i, genotype in enumerate( constants.genotypes ):
            genome_genotype = "_".join( ( genome, genotype ) )
            
            priors[genome]['location']['alpha'][i] = parser.getfloat( 'location_alpha', genome_genotype )
            priors[genome]['location']['beta'][i] = parser.getfloat( 'location_beta', genome_genotype )

            priors[genome]['precision']['shape'][i] = parser.getfloat( 'precision_shape', genome_genotype )
            priors[genome]['precision']['scale'][i] = parser.getfloat( 'precision_scale', genome_genotype )
            
    return priors

def get_independent_binomial_density_priors( parser, priors ):
    for genome in constants.genomes:
        priors[genome]['alpha'] = np.zeros( ( 3, ) )
        priors[genome]['beta'] = np.zeros( ( 3, ) ) 
        
        for i, genotype in enumerate( constants.genotypes ):
            genome_genotype = "_".join( ( genome, genotype ) )
            
            priors[genome]['alpha'][i] = parser.getfloat( 'alpha', genome_genotype )
            priors[genome]['beta'][i] = parser.getfloat( 'beta', genome_genotype )
        
    return priors
    
def parse_independent_parameters_file( parameters_file_name, density ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( parameters_file_name )
    
    parameters = {}
    parameters['pi'] = np.zeros( ( 2, 3 ) )
    
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        parameters['pi'][0, i] = parser.getfloat( 'pi', normal_genotype )
        parameters['pi'][1, i] = parser.getfloat( 'pi', tumour_genotype )
    
    parameters['pi'][0, :] = parameters['pi'][0, :] / parameters['pi'][0, :].sum()
    parameters['pi'][1, :] = parameters['pi'][1, :] / parameters['pi'][1, :].sum()
    
    if density == "beta_binomial":
        normal_parameters, tumour_parameters = get_independent_beta_binomial_density_parameters( parser, parameters )
    elif density == "binomial":
        normal_parameters, tumour_parameters = get_independent_binomial_density_parameters( parser, parameters )

    return normal_parameters, tumour_parameters

def get_independent_beta_binomial_density_parameters( parser, parameters ):
    parameters['alpha'] = np.zeros( ( 2, 3 ) )
    parameters['beta'] = np.zeros( ( 2, 3 ) )
    
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

def get_independent_binomial_density_parameters( parser, parameters ):
    parameters['mu'] = np.zeros( ( 2, 3 ) )
    
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        parameters['mu'][0, i] = parser.getfloat( 'mu', normal_genotype )
        parameters['mu'][1, i] = parser.getfloat( 'mu', tumour_genotype )
        
    normal_parameters = {'pi' : parameters[0, :], 'mu' : parameters[0, :]}
    tumour_parameters = {'pi' : parameters[1, :], 'mu' : parameters[1, :]}
        
    return normal_parameters, tumour_parameters

#=======================================================================================================================
# Joint Models
#=======================================================================================================================       

def write_chr_responsibilities( model, reader, writer, parameters, chr_name ):
    counts = reader.get_counts( chr_name )
    
    jcnt_rows = reader.get_rows( chr_name )
    
    end = reader.get_chr_size( chr_name )
    
    n = int( 1e5 )
    
    start = 0
    stop = min( n, end )
    
    while start < end:
        sub_counts = counts[start:stop]
        sub_rows = jcnt_rows[start:stop]
        
        data = JointData( sub_counts )
        responsibilities = model.classify( data, parameters )
        
        writer.write_data( chr_name, sub_rows, responsibilities )
        
        start = stop
        stop = min( stop + n, end )

def run_joint_model( args ):
    if args.density == "beta_binomial":
        model = JointBetaBinomialModel()
    elif args.density == "binomial":
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
        
        priors = parse_joint_priors_file( args.priors_file, data.nrows, args.density )
        
        writer.write_priors( priors )
        
        parameters = model.train( data, priors, args.max_iters, args.convergence_threshold )
    else:
        parameters = parse_joint_parameters_file( args.params_file, args.density )
        print parameters
    
    writer.write_parameters( parameters )
    
    chr_list = reader.get_chr_list()
    
    for chr_name in sorted( chr_list ):
        write_chr_responsibilities( model, reader, writer, parameters, chr_name )
    
    reader.close()
    writer.close()
        
def parse_joint_priors_file( priors_file_name, n, density ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( priors_file_name )
    
    priors = {}
    priors['kappa'] = np.zeros( ( 9, ) )
    
    for i, genotype_tuple in enumerate( constants.joint_genotypes ):
        genotype = "_".join( genotype_tuple )
        
        priors['kappa'][i] = parser.getfloat( 'kappa', genotype )
        
    scaling = parser.getint( 'kappa', 'scaling' )
    
    priors = scale_priors( priors, n, scaling )
    
    if density == "beta_binomial":
        priors = get_joint_beta_binomial_density_priors( parser, priors )
    elif density == "binomial":
        priors = get_joint_binomial_density_priors( parser, priors )
        
    return priors

def get_joint_beta_binomial_density_priors( parser, priors ):
    priors['location'] = np.zeros( ( 2, 3, 2 ) )
    priors['precision'] = np.zeros( ( 2, 3, 3 ) )
          
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        priors['location'][0, i, 0] = parser.getfloat( 'location_alpha', normal_genotype )
        priors['location'][0, i, 1] = parser.getfloat( 'location_beta', normal_genotype )
        
        priors['location'][1, i, 0] = parser.getfloat( 'location_alpha', tumour_genotype )
        priors['location'][1, i, 1] = parser.getfloat( 'location_beta', tumour_genotype )
        
        priors['precision'][0, i, 0] = parser.getfloat( 'precision_shape', normal_genotype )
        priors['precision'][0, i, 1] = parser.getfloat( 'precision_scale', normal_genotype )
        priors['precision'][0, i, 2] = parser.getfloat( 'precision_min', normal_genotype )
        
        priors['precision'][1, i, 0] = parser.getfloat( 'precision_shape', tumour_genotype )
        priors['precision'][1, i, 1] = parser.getfloat( 'precision_scale', tumour_genotype )
        priors['precision'][1, i, 2] = parser.getfloat( 'precision_min', tumour_genotype )
        
    return priors

def get_joint_binomial_density_priors( parser, priors ):
    priors['alpha'] = np.zeros( ( 2, 3 ) )
    priors['beta'] = np.zeros( ( 2, 3 ) )
    
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        priors['alpha'][0, i] = parser.getfloat( 'alpha', normal_genotype )
        priors['beta'][0, i] = parser.getfloat( 'beta', normal_genotype )
        
        priors['alpha'][1, i] = parser.getfloat( 'alpha', tumour_genotype )
        priors['beta'][1, i] = parser.getfloat( 'beta', tumour_genotype )
        
    return priors
    
def parse_joint_parameters_file( parameters_file_name, density ):
    parser = ConfigParser.SafeConfigParser()
    parser.read( parameters_file_name )
    
    parameters = {}
    parameters['pi'] = np.zeros( ( 9, ) )
    
    for i, genotype_tuple in enumerate( constants.joint_genotypes ):
        genotype = "_".join( genotype_tuple )
        
        parameters['pi'][i] = parser.getfloat( 'pi', genotype )
    
    parameters['pi'] = parameters['pi'] / parameters['pi'].sum()
    
    if density == "beta_binomial":
        get_joint_beta_binomial_density_paramters( parser, parameters )
    elif density == "binomial":
        get_joint_binomial_density_paramters( parser, parameters )
    
    return parameters


def get_joint_beta_binomial_density_paramters( parser, parameters ):
    parameters['alpha'] = np.zeros( ( 2, 3 ) )
    parameters['beta'] = np.zeros( ( 2, 3 ) )
    
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        parameters['alpha'][0, i] = parser.getfloat( 'alpha', normal_genotype )
        parameters['alpha'][1, i] = parser.getfloat( 'alpha', tumour_genotype )
        
        parameters['beta'][0, i] = parser.getfloat( 'beta', normal_genotype )
        parameters['beta'][1, i] = parser.getfloat( 'beta', tumour_genotype )
        
    return parameters

def get_joint_binomial_density_paramters( parser, parameters ):
    parameters['mu'] = np.zeros( ( 2, 3 ) )
            
    for i, genotype in enumerate( constants.genotypes ):
        normal_genotype = "_".join( ( 'normal', genotype ) )
        tumour_genotype = "_".join( ( 'tumour', genotype ) )
        
        parameters['mu'][0, i] = parser.getfloat( 'mu', normal_genotype )
        parameters['mu'][1, i] = parser.getfloat( 'mu', tumour_genotype )
        
    return parameters

#=======================================================================================================================
# Chromosome Fit
#=======================================================================================================================
def run_chromosome_model( args ):
    if args.density == "beta_binomial":
        model = JointBetaBinomialModel()
    elif args.density == "binomial":
        model = JointBinomialModel()
    
    reader = JointCountsReader( args.jcnt_file_name )
    writer = JointSnvMixWriter( args.jsm_file_name )
    
    n = 1
    
    priors = parse_joint_priors_file( args.priors_file, n, args.density )
        
    writer.write_priors( priors )
    
    chr_list = reader.get_chr_list()
    
    for chr_name in sorted( chr_list ):
        print chr_name
        
        if args.subsample_size > 0:
            counts = chrom_subsample( reader, chr_name, args.subsample_size )
        else:        
            counts = reader.get_counts( chr_name )
        
        data = JointData( counts )
        print data.nrows
        
        parameters = model.train( data, priors, args.max_iters, args.convergence_threshold )
        
        writer.write_chr_parameters( parameters, chr_name )
        
        print "Converged parameter for chromosome {0} are.".format( chr_name )
        print parameters
        
        write_chr_responsibilities( model, reader, writer, parameters, chr_name )
    
    reader.close()
    writer.close()

def chrom_subsample( reader, chr_name, sample_size ):
    chr_size = reader.get_chr_size( chr_name=chr_name )
    
    sample_size = min( chr_size, sample_size )
    
    chr_sample_indices = random.sample( xrange( chr_size ), sample_size )
    
    chr_counts = reader.get_counts( chr_name )
    
    chr_sample = chr_counts[chr_sample_indices]
    
    return chr_sample
