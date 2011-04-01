'''
Created on 2011-02-11

@author: Andrew Roth
'''
import math
import random
import multiprocessing

import numpy as np
from scipy.cluster.vq import kmeans2


from joint_snv_mix import constants
from joint_snv_mix.classification.data import JointData
from joint_snv_mix.classification.latent_variables import EMLatentVariables
from joint_snv_mix.classification.likelihoods import joint_beta_binomial_log_likelihood, joint_binomial_log_likelihood
from joint_snv_mix.classification.lower_bounds import EMLowerBound
from joint_snv_mix.classification.model_runners import ModelRunner
from joint_snv_mix.classification.models import EMModel, EMModelTrainer
from joint_snv_mix.classification.posteriors import EMPosterior
from joint_snv_mix.classification.utils.beta_binomial_map_estimators import get_mle_p
from joint_snv_mix.classification.utils.log_pdf import log_translated_gamma_pdf, log_beta_pdf
from joint_snv_mix.file_formats.cncnt import ConanCountsReader
from joint_snv_mix.file_formats.cnsm import ConanSnvMixWriter

def run_conan( args ):
    if args.density == 'binomial':
        runner = ConanBinomialRunner()
    elif args.density == 'beta_binomial':       
        runner = ConanBetaBinomialRunner()
    
    args.train = True
    
    runner.run( args )

#=======================================================================================================================
# Runner
#=======================================================================================================================
class ConanModelRunner( ModelRunner ):
    def __init__( self ):
        ModelRunner.__init__( self )
        
        self.data_class = JointData
        self.parameters = {}
        self.priors = {}
    
    def run( self, args ):
        self.reader = ConanCountsReader( args.cncnt_file_name )
        self.writer = ConanSnvMixWriter( args.cnsm_file_name )
        
        ModelRunner.run( self, args )
    
    def _classify( self, args ):
        cn_states = self.reader.get_cn_states()
        
        for cn_state in sorted( cn_states ):
            self._classify_cn_state( cn_state )
                
    def _classify_cn_state( self, cn_state ):        
        chr_list = self.reader.get_table_list( cn_state )

        for chr_name in sorted( chr_list ):            
            self._classify_chromosome( cn_state, chr_name )    
    
    def _train( self, args ):               
        cn_states = self.reader.get_cn_states()
                
        for cn_state in sorted( cn_states ):             
            self._train_cn_state( cn_state, args )
        
        self._write_priors()
    
    def _train_cn_state( self, cn_state, args ):                   
        if args.subsample_size > 0:
            counts = self._subsample( cn_state, args.subsample_size )
        else:
            counts = self.reader.get_counts()

        nclass = {}
        nclass['normal'] = 3
        nclass['tumour'] = constants.cn_state_map[cn_state]
        
        priors = self._get_priors( nclass )
                        
        model = self.model_class( nclass )
        
        data = self.data_class( counts )
        
        self.parameters[cn_state] = model.train( 
                                                data, priors,
                                                args.max_iters,
                                                args.convergence_threshold
                                                )
        
        self.priors[cn_state] = priors
    
    def _classify_chromosome( self, cn_state, chr_name ):
        nclass = {}
        nclass['normal'] = 3
        nclass['tumour'] = constants.cn_state_map[cn_state]              
        
        model = self.model_class( nclass )
                
        counts = self.reader.get_counts( cn_state, chr_name )
       
        jcnt_table = self.reader.get_table( cn_state, chr_name )
        
        end = self.reader.get_number_of_table_rows( cn_state, chr_name )

        n = int( 1e5 )
        start = 0
        stop = min( n, end )
        

        while start < end:
            sub_counts = counts[start:stop]
            sub_rows = jcnt_table[start:stop]
                              
            data = self.data_class( sub_counts )            
                
            resp = model.classify( data, self.parameters[cn_state] )
        
            self.writer.write_data( cn_state, chr_name, sub_rows, resp )
            
            start = stop
            stop = min( stop + n, end )
            
    def _subsample( self, cn_state, sample_size ):
        chr_list = self.reader.get_table_list( cn_state )
        
        sample = []
        
        nrows = self.reader.get_data_set_size( cn_state )
        
        for chr_name in chr_list:
            chr_size = self.reader.get_number_of_table_rows( cn_state, chr_name )
            
            chr_sample_size = math.floor( float( chr_size ) / nrows * sample_size )
            
            chr_sample_size = int( chr_sample_size )
            
            chr_sample_size = min( chr_size, chr_sample_size )
            
            chr_sample_indices = random.sample( xrange( chr_size ), chr_sample_size )
            
            chr_counts = self.reader.get_counts( cn_state, chr_name )
            
            chr_sample = chr_counts[chr_sample_indices]
            
            sample.append( chr_sample )
            
        sample = np.vstack( sample )
        
        return sample
    
class ConanBetaBinomialRunner( ConanModelRunner ):
    def __init__( self ):
        ConanModelRunner.__init__( self )
        
        self.model_class = ConanBetaBinomialModel
    
    def _get_priors( self, nclass ):
        ncomponents = nclass['normal'] * nclass['tumour']
        
        priors = {}                
        priors['kappa'] = 10 * np.ones( ( ncomponents, ) )
        
        for genome in constants.genomes:
            priors[genome] = {}
            
            priors[genome]['precision'] = {}
            
            priors[genome]['precision']['shape'] = 20 * np.ones( ( nclass[genome], ) )
            priors[genome]['precision']['scale'] = 50 * np.ones( ( nclass[genome], ) )
            priors[genome]['precision']['min'] = 2 * np.ones( ( nclass[genome], ) )
            
            priors[genome]['location'] = {}
            priors[genome]['location']['alpha'] = np.linspace( 1000, 1, nclass[genome] )
            priors[genome]['location']['beta'] = np.linspace( 1, 1000, nclass[genome] )
        
        return priors

class ConanBinomialRunner( ConanModelRunner ):
    def __init__( self ):
        ConanModelRunner.__init__( self )
        
        self.model_class = ConanBinomialModel
    
    def _get_priors( self, nclass ):
        ncomponents = nclass['normal'] * nclass['tumour']
        
        priors = {}                
        priors['kappa'] = 10 * np.ones( ( ncomponents, ) )
        
        for genome in constants.genomes:
            priors[genome] = {}
            priors[genome]['mu'] = {}
            
            priors[genome]['mu']['alpha'] = np.linspace( 1000, 2, nclass[genome] )
            priors[genome]['mu']['beta'] = np.linspace( 2, 1000, nclass[genome] )
        
        return priors            

#=======================================================================================================================
# Model
#=======================================================================================================================
class ConanBetaBinomialModel( EMModel ):
    def __init__( self, nclass ):
        self.trainer_class = ConanBetaBinomialModelTrainer
        
        self.log_likelihood_func = joint_beta_binomial_log_likelihood
        
        self.nclass = nclass
        
    def train( self, data, priors, max_iters, tolerance ):
        '''
        Train the model using EM.
        
        Input: JointData object
        '''   
        trainer = self.trainer_class( data, self.nclass, max_iters, tolerance, priors )
        
        parameters = trainer.run()
        
        trainer.responsibilities = []
                
        return parameters
    
class ConanBinomialModel( EMModel ):
    def __init__( self, nclass ):
        self.trainer_class = ConanBinomialModelTrainer
        
        self.log_likelihood_func = joint_binomial_log_likelihood
        
        self.nclass = nclass
        
    def train( self, data, priors, max_iters, tolerance ):
        '''
        Train the model using EM.
        
        Input: JointData object
        '''   
        trainer = self.trainer_class( data, self.nclass, max_iters, tolerance, priors )
        
        parameters = trainer.run()
        
        trainer.responsibilities = []
                
        return parameters

#=======================================================================================================================
# Model Trainers
#=======================================================================================================================
class ConanBetaBinomialModelTrainer( EMModelTrainer ):
    def __init__( self, data, nclass, max_iters, tolerance, priors ):
        self.nclass = nclass
        
        EMModelTrainer.__init__( self, data, max_iters, tolerance, priors )        
        
    def _init_components( self ):
        self.latent_variables = ConanBetaBinomialLatentVariables( self.data, self.nclass )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = ConanBetaBinomialPosterior( self.data, self.priors,
                                                     self.responsibilities, self.nclass )
        
        self.lower_bound = ConanBetaBinomialLowerBound( self.data, self.priors )
        
class ConanBinomialModelTrainer( EMModelTrainer ):
    def __init__( self, data, nclass, max_iters, tolerance, priors ):
        self.nclass = nclass
        
        EMModelTrainer.__init__( self, data, max_iters, tolerance, priors )        
        
    def _init_components( self ):
        self.latent_variables = ConanBinomialLatentVariables( self.data, self.nclass )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = ConanBinomialPosterior( self.data, self.priors,
                                                 self.responsibilities, self.nclass )
        
        self.lower_bound = ConanBinomialLowerBound( self.data, self.priors )
        
#=======================================================================================================================
# Latent Variables
#=======================================================================================================================
class ConanLatentVariables( EMLatentVariables ):
    def __init__( self, data, nclass ):
        self.nclass = nclass
        self.ncomponents = self.nclass['normal'] * self.nclass['tumour']
        
        EMLatentVariables.__init__( self, data )
    
    def _init_responsibilities( self, data ):
        '''
        Intialise responsibilities via k-means clustering.
        '''
        shape = ( data.nrows, self.ncomponents )
        
        responsibilities = np.zeros( shape )
        
        labels = {}
        for genome in constants.genomes:
            a = np.asarray( data.a[genome], dtype=np.float64 )
            b = np.asarray( data.b[genome], dtype=np.float64 )
            d = a + b
            p = a / d
              
            init_centers = np.linspace( 1, 0, self.nclass[genome] )                    
            
            clustering_result = kmeans2( p, init_centers, minit='matrix' )
            
            labels[genome] = clustering_result[1]
            
            print "Initial class ceneters : ", clustering_result[0]

        labels = self.nclass['normal'] * labels['normal'] + labels['tumour']

        for id in range( self.ncomponents ):
            indices = ( labels == id )
            
            responsibilities[indices, id] = 1.
        
        self.responsibilities = responsibilities
        
class ConanBetaBinomialLatentVariables( ConanLatentVariables ):
    def __init__( self, data, nclass ):
        ConanLatentVariables.__init__( self, data, nclass )
        
        self.likelihood_func = joint_beta_binomial_log_likelihood
        
class ConanBinomialLatentVariables( ConanLatentVariables ):
    def __init__( self, data, nclass ):
        ConanLatentVariables.__init__( self, data, nclass )
        
        self.likelihood_func = joint_binomial_log_likelihood

#=======================================================================================================================
# Posteriors
#=======================================================================================================================
class ConanPosterior( EMPosterior ):
    def _get_marginals( self ):        
        nclass = self.nclass
        nrows = self.data.nrows
        responsibilities = self.responsibilities
    
        shape = ( nrows, nclass['normal'], nclass['tumour'] )
        
        responsibilities = responsibilities.reshape( shape )
        
        marginals = {}
        
        marginals['normal'] = responsibilities.sum( axis=2 )
        marginals['tumour'] = responsibilities.sum( axis=1 )
    
        return marginals

class ConanBetaBinomialPosterior( ConanPosterior ):
    def __init__( self, data, priors, responsibilities, nclass ):
        self.nclass = nclass
        self.ncomponents = self.nclass['normal'] * self.nclass['tumour']
               
        EMPosterior.__init__( self, data, priors, responsibilities )
    
    def _init_parameters( self ):
        '''
        Initialise parameters. This is only necessary to initialise gradient descent. 
        '''
        self.parameters = {}

        self._update_mix_weights()

        for genome in constants.genomes:
            self.parameters[genome] = {}
            
            self.parameters[genome]['alpha'] = np.array( [100] * self.nclass[genome],
                                                 dtype=np.float64 )
            
            self.parameters[genome]['beta'] = np.array( [100] * self.nclass[genome],
                                                 dtype=np.float64 )


    
    def _update_density_parameters( self ):        
        marginals = self._get_marginals()
        
        print "Begining numerical optimisation of alpha and beta."
        
        vars = []
        
        for genome in constants.genomes:
            for component in range( self.nclass[genome] ):
                a = self.data.a[genome]
                b = self.data.b[genome]
                
                x = np.zeros( ( 2, ) )
                
                x[0] = self.parameters[genome]['alpha'][component]
                x[1] = self.parameters[genome]['beta'][component]
                
                resp = marginals[genome][:, component]
                
                precision_prior = self.priors[genome]['precision']
                location_prior = self.priors[genome]['location']
                
                vars.append( [x, a, b, resp, location_prior, precision_prior, component] )

        pool = multiprocessing.Pool( processes=self.ncomponents, maxtasksperchild=1 )
        results = pool.map( get_mle_p, vars )
        pool.close()
        
        i = 0
        for genome in constants.genomes:
            for component in range( self.nclass[genome] ):                
                self.parameters[genome]['alpha'][component] = results[i][0]
                self.parameters[genome]['beta'][component] = results[i][1]
                
                i += 1

class ConanBinomialPosterior( ConanPosterior ):
    def __init__( self, data, priors, responsibilities, nclass=3 ):
        self.nclass = nclass
        
        EMPosterior.__init__( self, data, priors, responsibilities )
    
    def _init_parameters( self ):
        self.parameters = {}
        self.parameters['normal'] = {}
        self.parameters['tumour'] = {}
        
        self._update_mix_weights()
        
        self._update_density_parameters()
           
    def _update_density_parameters( self ):
        marginals = self._get_marginals()
        
        for genome in constants.genomes:
            a = self.data.a[genome]
            b = self.data.b[genome]
            
            alpha = self.priors[genome]['mu']['alpha']
            beta = self.priors[genome]['mu']['beta']
            
            tau = marginals[genome]
            
            self.parameters[genome]['mu'] = self._update_mu( a, b, alpha, beta, tau )
    
    def _update_mu( self, a, b, alpha, beta, tau ):       
        d = a + b
        
        n = self.data.nrows
        shape = ( n, 1 )
        
        a = a.reshape( shape )
        d = d.reshape( shape )
        
        ref_sum = np.sum( tau * a, axis=0 )
        
        depth_sum = np.sum( tau * d, axis=0 )
        
        numerator = ref_sum + alpha - 1.
        
        denominator = depth_sum + alpha + beta - 2.
        
        return np.exp( np.log( numerator ) - np.log( denominator ) )

#===================================================================================================================
# Lower Bounds
#===================================================================================================================
class ConanBetaBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_beta_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        precision_term = 0
        location_term = 0
                
        for genome in constants.genomes:   
            alpha = self.parameters[genome]['alpha']
            beta = self.parameters[genome]['beta']
            
            s = alpha + beta
            mu = alpha / s
    
            precision_priors = self.priors[genome]['precision']
            location_priors = self.priors[genome]['location']
            
            precision_term += np.sum( log_translated_gamma_pdf( s,
                                                                precision_priors['shape'],
                                                                precision_priors['scale'],
                                                                precision_priors['min'] 
                                                                ) )
            
            location_term += np.sum( log_beta_pdf( mu,
                                                   location_priors['alpha'],
                                                   location_priors['beta'] ) )
                
        return precision_term + location_term

class ConanBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        log_prior = 0.
        
        for genome in constants.genomes:
            mu = self.parameters[genome]['mu']
            
            alpha = self.priors[genome]['mu']['alpha']
            beta = self.priors[genome]['mu']['beta']
            
            log_prior += np.sum( log_beta_pdf( mu, alpha, beta ) )

        return log_prior
