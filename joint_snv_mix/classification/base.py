'''
Abstract classes which all models for classifying paired data should sub-class. 

Created on 2011-03-31

@author: Andrew Roth
'''
from ConfigParser import ConfigParser
import math

import numpy as np

from joint_snv_mix import constants

from joint_snv_mix.file_formats.jcnt import JointCountsReader
from joint_snv_mix.file_formats.jsm import JointSnvMixWriter

from joint_snv_mix.classification.utils.normalise import log_space_normalise_rows
from joint_snv_mix.classification.utils.log_pdf import log_dirichlet_pdf
from joint_snv_mix.classification.utils.data import JointData

class ProbabilisticModelRunner(object):
    '''
    Class for running a probabilistic model for paired data analysis.
    '''
    def __init__(self, model, priors_parser, parameters_parser):
        self.model = model
        self.priors_parser = priors_parser
        self.parameter_parser = parameters_parser
    
    def run(self, args):
        '''
        Run a full analysis from arguments.
        '''
        self.reader = JointCountsReader(args.jcnt_file_name)
        self.writer = JointSnvMixWriter(args.jsm_file_name)
          
        # Load parameters by training or from file.
        if args.train:
            self._train(args)
        else:
            self._load_parameters(args)
        
        self._write_parameters()
        
        self._classify(args)
                            
        self.reader.close()
        self.writer.close()
    
    def _train(self, args):        
        if args.subsample_size > 0:
            counts = self._subsample(
                                     args.subsample_size, 
                                     args.min_train_depth, 
                                     args.max_train_depth
                                     )
        else:
            counts = self.reader.get_counts()
                   
        self.priors_parser.load_from_file(args.priors_file)
        self.priors = self.priors_parser.to_dict()
        
        self._write_priors()
        
        data = JointData(counts)
        
        self.parameters = self.model.train(data,
                                           self.priors,
                                           args.max_iters,
                                           args.convergence_threshold)
    
    def _classify(self, args):
        chr_list = self.reader.get_table_list()
        
        for chr_name in sorted(chr_list):
            self._classify_chromosome(chr_name)
            
    def _classify_chromosome(self, chrom):
        counts = self.reader.get_counts(chrom)
        jcnt_table = self.reader.get_table(chrom)
        
        end = self.reader.get_number_of_table_rows(chrom)

        n = int(1e6)
        start = 0
        stop = min(n, end)
        
        # Classify using blocking for speedup.
        while start < end:
            sub_counts = counts[start:stop]
            sub_jcnt_rows = jcnt_table[start:stop]
                              
            data = JointData(sub_counts)            
                
            resp = self.model.classify(data, self.parameters)
            
            sub_jcnt_rows = sub_jcnt_rows.tolist()
            resp = resp.tolist()

            for jcnt_row, probs in zip(sub_jcnt_rows, resp):
                jsm_row = []
                jsm_row.extend(jcnt_row)
                jsm_row.extend(probs)
                
                self.writer.add_row(chrom, jsm_row)
            
            print "Classifying row {0} out of {1} for chromosome {2}".format(stop, end, chrom)
            print "\t".join([str(x) for x in jsm_row])
            
            start = stop
            stop = min(stop + n, end)
            
    def _load_parameters(self, args):
        self.parameter_parser.load_from_file(args.params_file)
        
        self.parameters = self.parameter_parser.to_dict()
            
    def _write_parameters(self):
        self.writer.write_parameters(self.parameters)
        
    def _write_priors(self):
        self.writer.write_priors(self.priors)
            
    def _subsample(self, sample_size, min_depth, max_depth):
        table_list = self.reader.get_table_list()
        
        sample = []
        
        nrows = self.reader.get_data_set_size()
        
        for chrom in table_list:
            table_nrows = self.reader.get_number_of_table_rows(chrom=chrom)
            
            chrom_sample_size = math.ceil(float(table_nrows) / nrows * sample_size)
            
            chrom_sample_size = int(chrom_sample_size)
            
            chrom_sample = self.reader.get_random_counts_subsample(
                                                                   chrom, 
                                                                   chrom_sample_size, 
                                                                   min_depth=min_depth,
                                                                   max_depth=max_depth
                                                                   )
            
            print np.min(chrom_sample.sum(axis=1))
            print np.max(chrom_sample.sum(axis=1))
            
            sample.append(chrom_sample)
            
        sample = np.vstack(sample)
        
        return sample
    
#=======================================================================================================================
# EM Related Classes
#=======================================================================================================================
class EMModel(object):
    '''
    Abstract class for probabilistic models fit via EM.
    '''
    def __init__(self):
        self.trainer_class = None
        self.log_likelihood_func = None
            
    def train(self, data, priors, max_iters, tolerance):
        '''
        Train the model using EM to find the MAP estimate for parameters.
        '''   
        trainer = self.trainer_class(data, max_iters, tolerance, priors)
        
        parameters = trainer.run()
        
        trainer.responsibilities = []
                
        return parameters

    def classify(self, data, parameters):
        '''
        Classify the given data based on the parameters passed in.
        '''        
        log_responsibilities = self.log_likelihood_func(data, parameters)
        
        responsibilities = log_space_normalise_rows(log_responsibilities)
        
        return responsibilities

class EMModelTrainer(object):
    '''
    Abstract class for performing the actual training of EM model.
    '''
    def __init__(self, data, max_iters, tolerance, priors):
        self.max_iters = max_iters
        
        self.tolerance = tolerance
        
        self.data = data

        self.priors = priors
            
        self._init_components()
        
    def run(self):
        iters = 0
        converged = False
        
        parameters = self.posterior.parameters
        old_posterior_value = self.lower_bound.get_lower_bound(parameters)
  
        while not converged:
            self._M_step()
            self._E_step()

            posterior_value = self.lower_bound.get_lower_bound(self.parameters)
            
            if iters > 0:
                posterior_change = (posterior_value - old_posterior_value) / abs(old_posterior_value)
            else:
                posterior_change = float('inf')
            
            self._print_diagnostic_message(iters, posterior_value, old_posterior_value, posterior_change)
            old_posterior_value = posterior_value
            
            if posterior_change < 0:
                print "Posterior decreased. This could be a bug or overly stringent convergence criterion."
                converged = True
            
            elif posterior_change < self.tolerance:
                converged = True
                            
            if iters >= self.max_iters:
                print "Maximum numbers of EM iterations exceeded. Exiting training."                
                converged = True
            
            iters += 1
                     
        return self.parameters
                  
    def _E_step(self):
        self.latent_variables.update(self.parameters)
        self.responsibilities = self.latent_variables.responsibilities
        
    def _M_step(self):
        self.posterior.update(self.responsibilities)
        self.parameters = self.posterior.parameters
        
    def _print_diagnostic_message(self, iters, posterior_value, old_posterior_value, posterior_change):
        print "#" * 100
        print "# Diagnostics."
        print "#" * 100
        print "Number of iterations : ", iters
        print "New posterior : ", posterior_value
        print "Old posterior : ", old_posterior_value 
        print "Posterior change : ", posterior_change
    
        print "Parameters :"
        
        for param_name, param_value in self.posterior.parameters.items():
            print param_name, param_value     
    
    def _init_components(self):
        raise NotImplemented

class EMLatentVariables(object):
    '''
    Abstract class for storing and updating the expected value of the latent variables (responsibilities) in EM models.
    '''    
    def __init__(self, data):       
        self.data = data
        
        self._init_responsibilities(data)
        
        self.likelihood_func = None

    def update(self, parameters):       
        log_responsibilities = self.likelihood_func(self.data, parameters)

        self.responsibilities = log_space_normalise_rows(log_responsibilities)
       
    def _init_responsibilities(self, data):
        NotImplementedError     
           
class EMPosterior(object):
    '''
    Class for updating the EM model parameters.
    '''
    def __init__(self, data, priors, responsibilities):
        self.data = data
        self.priors = priors
        self.responsibilities = responsibilities
        
        self._init_parameters()
    
    def _init_parameters(self):
        raise NotImplemented

    def update(self, responsibilities):
        self.responsibilities = responsibilities     

        self._update_mix_weights()
        
        self._update_density_parameters()

    def _update_mix_weights(self):
        N_g = self.responsibilities.sum(axis=0)

        mix_weights = N_g + self.priors['kappa'] - 1

        mix_weights = np.exp(np.log(mix_weights) - np.log(mix_weights.sum()))

        self.parameters['pi'] = mix_weights
        
    def _update_density_parameters(self):
        raise NotImplemented

class EMLowerBound(object):
    '''
    Class for computing the lower bound (log posterior) for an EM model.
    '''
    def __init__(self, data, priors):                
        self.priors = priors
        self.data = data
        
        self.log_likelihood_func = None

    def get_lower_bound(self, parameters):
        self.parameters = parameters
        
        log_likelihood = self._get_log_likelihood()

        log_mix_weight_prior = self._get_log_mix_weight_prior()

        log_density_parameters_prior = self._get_log_density_parameters_prior()

        lower_bound = log_likelihood + log_mix_weight_prior + log_density_parameters_prior

        return lower_bound
    
    def _get_log_likelihood(self):
        log_likelihoods = self.log_likelihood_func(self.data, self.parameters)
        
        log_likelihood = np.logaddexp.reduce(log_likelihoods, axis=1).sum()

        return log_likelihood

    def _get_log_mix_weight_prior(self):
        pi = self.parameters['pi']
        kappa = self.priors['kappa']
        
        log_prior = log_dirichlet_pdf(pi, kappa)
        log_prior = log_prior.sum()

        return log_prior

    def _get_log_density_parameters_prior(self):
        raise NotImplemented

class PriorParser(object):
    '''
    Abstract class for loading hyper-parameters in prior distribution. This is only relevant for EM based models that
    are going to be trained.
    '''    
    def __init__(self):
        self.nclass = {}
        self.nclass['normal'] = 3
        self.nclass['tumour'] = 3        
                
        self.priors = {}
        
        for genome in constants.genomes:
            self.priors[genome] = {}
    
    def load_from_file(self, file_name):                
        self.parser = ConfigParser()
        self.parser.read(file_name)
        
        self._load_mix_weight_priors()
        self._load_density_priors()
        
    def to_dict(self):
        return self.priors
        
    def _load_mix_weight_priors(self):       
        raise NotImplemented
    
    def _load_density_priors(self):        
        for genome in constants.genomes:
            for param_name in self.parameter_names:
                self.priors[genome][param_name] = {}
                    
                for hyper_param_name in self.hyper_parameter_names[param_name]:
                    self.priors[genome][param_name][hyper_param_name] = np.zeros((self.nclass[genome],))
                    
                    self._load_hyperparameter(genome, param_name, hyper_param_name)
                    
    def _load_hyperparameter(self, genome, param_name, hyper_param_name):                           
        for i, genotype in enumerate(constants.genotypes):                
            genome_genotype = "_".join((genome, genotype))
            
            section_name = "_".join((param_name, hyper_param_name))
            
            self.priors[genome][param_name][hyper_param_name][i] = self.parser.getfloat(section_name, genome_genotype)
    
class ParameterParser(object):
    '''
    Abstract class for loading model parameters from file. This is only relevant to EM based models that are not
    training parameters.
    ''' 
    def __init__(self):
        self.nclass = {}
        self.nclass['normal'] = 3 
        self.nclass['tumour'] = 3        
                
        self.parameters = {}
        
        for genome in constants.genomes:
            self.parameters[genome] = {}
    
    def load_from_file(self, file_name):                
        self.parser = ConfigParser()
        self.parser.read(file_name)
        
        self._load_mix_weights()
        self._load_density_parameters()
        
    def to_dict(self):
        return self.parameters
        
    def _load_mix_weights(self):       
        raise NotImplemented
    
    def _load_density_parameters(self):        
        for genome in constants.genomes:
            for param_name in self.parameter_names:
                self.parameters[genome][param_name] = np.zeros((self.nclass[genome],))
                
                self._load_parameter(genome, param_name)
                    
    def _load_parameter(self, genome, param_name):                           
        for i, genotype in enumerate(constants.genotypes):                
            genome_genotype = "_".join((genome, genotype))
            
            self.parameters[genome][param_name][i] = self.parser.getfloat(param_name, genome_genotype)
