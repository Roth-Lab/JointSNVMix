'''
Created on 2011-02-03

@author: Andrew Roth

Classes for parsing parameters from config files.
'''
from ConfigParser import ConfigParser

import numpy as np

from joint_snv_mix import constants

class ParameterParser( object ):
    def __init__( self ):
        self.nclass = {}
        self.nclass['normal'] = 3
        self.nclass['tumour'] = 3        
                
        self.parameters = {}
        
        for genome in constants.genomes:
            self.parameters[genome] = {}
    
    def load_from_file( self, file_name ):                
        self.parser = ConfigParser()
        self.parser.read( file_name )
        
        self._load_mix_weights()
        self._load_density_parameters()
        
    def to_dict( self ):
        return self.parameters
        
    def _load_mix_weights( self ):       
        raise NotImplemented
    
    def _load_density_parameters( self ):        
        for genome in constants.genomes:
            for param_name in self.parameter_names:
                self.parameters[genome][param_name] = np.zeros( ( self.nclass[genome], ) )
                
                self._load_parameter( genome, param_name )
                    
    def _load_parameter( self, genome, param_name ):                           
        for i, genotype in enumerate( constants.genotypes ):                
            genome_genotype = "_".join( ( genome, genotype ) )
            
            self.parameters[genome][param_name][i] = self.parser.getfloat( param_name, genome_genotype )
            
#=======================================================================================================================
# Independent Models
#=======================================================================================================================
class IndepedendentParameterParser( ParameterParser ):    
    def _load_mix_weights( self ):       
        for genome in constants.genomes:            
            pi = np.zeros( ( self.nclass[genome], ) )
            
            for i, genotype in enumerate( constants.genotypes ):                
                genome_genotype = "_".join( ( genome, genotype ) )
            
                pi[i] = self.parser.getfloat( 'pi', genome_genotype )
                
            self.parameters[genome]['pi'] = pi / pi.sum()
                
class IndependentBinomialParameterParser( IndepedendentParameterParser ):
    def __init__( self ):
        IndepedendentParameterParser.__init__( self )
        
        self.parameter_names = ( 'mu', )

class IndependentBetaBinomialParameterParser( IndepedendentParameterParser ):
    def __init__( self ):
        IndepedendentParameterParser.__init__( self )
        
        self.parameter_names = ( 'alpha', 'beta' )
        
#=======================================================================================================================
# Joint Models
#=======================================================================================================================
class JointParameterParser( ParameterParser ):
    def __init__( self ):
        ParameterParser.__init__( self )
        
        self.ncomponent = sum( [self.nclass[genome] for genome in constants.genomes] )

    def _load_mix_weights( self ):       
        pi = np.zeros( ( self.ncomponent, ) )
            
        for i, genotype_tuple in enumerate( constants.joint_genotypes ):
            genotype = "_".join( genotype_tuple )
        
            pi[i] = self.parser.getfloat( 'pi', genotype )
            
        pi[9] = self.parser.getfloat( 'pi', 'junk' )
            
        self.parameters['pi'] = pi / pi.sum()
            
class JointBinomialParameterParser( JointParameterParser ):
    def __init__( self ):
        JointParameterParser.__init__( self )
        
        self.parameter_names = ( 'mu', )

class JointBetaBinomialParameterParser( JointParameterParser ):
    def __init__( self ):
        JointParameterParser.__init__( self )
        
        self.parameter_names = ( 'alpha', 'beta' )
        
if __name__ == "__main__":
    def print_params( file_name, parser ):
        parser.load_from_file( file_name )
    
        print parser.to_dict()

    bb_params_file = "../../config/joint_bb.params.cfg"    
    bb_parser = JointBetaBinomialParameterParser()
    print_params( bb_params_file, bb_parser )

    bin_params_file = "../../config/joint_bin.params.cfg"    
    bin_parser = JointBinomialParameterParser()
    print_params( bin_params_file, bin_parser )
    
    indep_bin_params_file = "../../config/indep_bin.params.cfg"    
    indep_bin_parser = IndependentBinomialParameterParser()
    print_params( indep_bin_params_file, indep_bin_parser )
    
    indep_bb_params_file = "../../config/indep_bb.params.cfg"    
    indep_bb_parser = IndependentBetaBinomialParameterParser()
    print_params( indep_bb_params_file, indep_bb_parser )
