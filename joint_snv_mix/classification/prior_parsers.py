'''
Created on 2011-02-03

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix import constants
from ConfigParser import ConfigParser

class PriorParser( object ):
    def __init__( self ):
        self.nclass = {}
        self.nclass['normal'] = 3
        self.nclass['tumour'] = 3       
        self.nclass['junk'] = 1 
                
        self.priors = {}
        
        for genome in constants.genomes:
            self.priors[genome] = {}
    
    def load_from_file( self, file_name ):                
        self.parser = ConfigParser()
        self.parser.read( file_name )
        
        self._load_mix_weight_priors()
        self._load_density_priors()
        
    def to_dict( self ):
        return self.priors
        
    def _load_mix_weight_priors( self ):       
        raise NotImplemented
    
    def _load_density_priors( self ):        
        for genome in constants.genomes:
            for param_name in self.parameter_names:
                self.priors[genome][param_name] = {}
                    
                for hyper_param_name in self.hyper_parameter_names[param_name]:
                    self.priors[genome][param_name][hyper_param_name] = np.zeros( ( self.nclass[genome], ) )
                    
                    self._load_hyperparameter( genome, param_name, hyper_param_name )
                    
    def _load_hyperparameter( self, genome, param_name, hyper_param_name ):
        if genome == "junk":
            
                                 
        for i, genotype in enumerate( constants.genotypes ):                
            genome_genotype = "_".join( ( genome, genotype ) )
            
            section_name = "_".join( ( param_name, hyper_param_name ) )
            
            self.priors[genome][param_name][hyper_param_name][i] = self.parser.getfloat( section_name, genome_genotype )

#=======================================================================================================================
# Independent Models
#=======================================================================================================================
class IndepedendentPriorParser( PriorParser ):    
    def _load_mix_weight_priors( self ):       
        for genome in constants.genomes:            
            self.priors[genome]['kappa'] = np.zeros( ( self.nclass[genome], ) )
            
            for i, genotype in enumerate( constants.genotypes ):                
                genome_genotype = "_".join( ( genome, genotype ) )
            
                self.priors[genome]['kappa'][i] = self.parser.getfloat( 'kappa', genome_genotype )
                
class IndependentBinomialPriorParser( IndepedendentPriorParser ):
    def __init__( self ):
        IndepedendentPriorParser.__init__( self )
        
        self.parameter_names = ( 'mu', )
        
        self.hyper_parameter_names = {}
        self.hyper_parameter_names['mu'] = ( 'alpha', 'beta' )

class IndependentBetaBinomialPriorParser( IndepedendentPriorParser ):
    def __init__( self ):
        IndepedendentPriorParser.__init__( self )
        
        self.parameter_names = ( 'location', 'precision' )
        
        self.hyper_parameter_names = {}
        self.hyper_parameter_names['location'] = ( 'alpha', 'beta' )
        self.hyper_parameter_names['precision'] = ( 'shape', 'scale', 'min' )
    
#=======================================================================================================================
# Joint Models
#=======================================================================================================================
class JointModelPriorParser( PriorParser ):
    def __init__( self ):
        PriorParser.__init__( self )
        
        self.ncomponent = sum( [self.nclass[genome] for genome in constants.genomes] )
    
    def _load_mix_weight_priors( self ):       
        self.priors['kappa'] = np.zeros( ( self.ncomponent, ) )
            
        for i, genotype_tuple in enumerate( constants.joint_genotypes ):
            if isinstance(genotype_tuple, tuple):
                genotype = "_".join( genotype_tuple )
            else:
                genotype = genotype_tuple
        
            self.priors['kappa'][i] = self.parser.getfloat( 'kappa', genotype )
        
class JointBinomialPriorParser( JointModelPriorParser ):
    def __init__( self ):
        JointModelPriorParser.__init__( self )
        
        self.parameter_names = ( 'mu', )
        
        self.hyper_parameter_names = {}
        self.hyper_parameter_names['mu'] = ( 'alpha', 'beta' )

class JointBetaBinomialPriorParser( JointModelPriorParser ):
    def __init__( self ):
        JointModelPriorParser.__init__( self )
        
        self.parameter_names = ( 'location', 'precision' )
        
        self.hyper_parameter_names = {}
        self.hyper_parameter_names['location'] = ( 'alpha', 'beta' )
        self.hyper_parameter_names['precision'] = ( 'shape', 'scale', 'min' )        

if __name__ == "__main__":
    def print_params( file_name, parser ):
        parser.load_from_file( file_name )
    
        print parser.to_dict()
    
    bb_file_name = "../../config/joint_bb.priors.cfg"
    bb_prior = JointBetaBinomialPriorParser()
    print_params( bb_file_name, bb_prior )
    
    bin_file_name = "../../config/joint_bin.priors.cfg"    
    bin_prior = JointBinomialPriorParser()    
    print_params( bin_file_name, bin_prior )
    
    indep_bin_file_name = "../../config/indep_bin.priors.cfg"
    indep_bin_prior = IndependentBinomialPriorParser()
    print_params( indep_bin_file_name, indep_bin_prior )

    indep_bb_file_name = "../../config/indep_bb.priors.cfg"
    indep_bb_prior = IndependentBetaBinomialPriorParser()
    print_params( indep_bb_file_name, indep_bb_prior )
