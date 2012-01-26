'''
Created on 2011-08-11

@author: Andrew Roth
'''
import sys

# Python imports
from joint_snv_mix.counter import JointBinaryCounter
from joint_snv_mix.samtools import BamFile, FastaFile

from joint_snv_mix.models.binomial import BinomialParameters, BinomialPriors, BinomialModel
from joint_snv_mix.models.beta_binomial import BetaBinomialParameters, BetaBinomialPriors, BetaBinomialModel
from joint_snv_mix.models.snv_mix_two import SnvMixTwoModel

# Cython imports
from joint_snv_mix.counter_row cimport JointBinaryCounterRow, JointBinaryData
from joint_snv_mix.models.abstract cimport MixtureModel
from joint_snv_mix.results cimport CResultsWriter
from joint_snv_mix.positions_counter cimport PositionsCounter

#=======================================================================================================================
# Functions for running classification.
#=======================================================================================================================
def classify(args):
    factory = ModelFactory()

    counter = factory.get_counter(args)
    model = factory.get_model(args) 
    
    classify_data_set(counter, model, args)

def classify_data_set(counter, MixtureModel model, args):
    cdef int tumour_var_counts
    cdef bint print_all
    cdef double somatic_threshold, p_somatic
    
    cdef JointBinaryCounterRow row
    cdef JointBinaryData data
    
    cdef CResultsWriter writer
    
    print_all = args.print_all_positions
    somatic_threshold = args.somatic_threshold
    
    if args.positions_file is not None:
        positions_counter = PositionsCounter(args.positions_file, counter)
    
    writer = CResultsWriter(file_name=args.out_file)
    
    if args.chromosome is None:
        if args.positions_file is None:
            refs = counter.refs
        else:
            refs = positions_counter.refs
    else:
        refs = [args.chromosome, ]
        
    for ref in sorted(refs):
        positions_iter = counter.get_ref_iterator(ref)
        
        if args.positions_file is not None:
            positions_iter = positions_counter.get_ref_iterator(positions_iter)
        
        for row in positions_iter:
            data = row._data
            tumour_var_counts = row._data._b_T
            
            if not print_all and tumour_var_counts == 0:
                continue
            
            model._predict(data)
            
            p_somatic = model._resp[1] + model._resp[2]
            
            if p_somatic < somatic_threshold:
                continue
            
            writer.write_position(row, model._resp)
    
    writer.close()
         
#=======================================================================================================================
# Functions for training a model.
#=======================================================================================================================
def train(args):
    factory = ModelFactory()

    counter = factory.get_counter(args)
    model = factory.get_model(args)     
    
    data_set = create_training_data_set(counter, args)
    
    model.fit(data_set, verbose=True)
    
    model.params.write_to_file(args.estimated_parameters_file)        

def create_training_data_set(counter, args):
    cdef int i, d_N, d_T, min_d_N, min_d_T, max_d_N, max_d_T, skip_size
    
    cdef JointBinaryCounterRow row
    cdef JointBinaryData data
    
    data_set = []
    
    i = 0
    
    min_d_N = args.min_normal_depth
    min_d_T = args.min_tumour_depth
    
    max_d_N = args.max_normal_depth
    max_d_T = args.max_tumour_depth
    
    skip_size = args.skip_size
    
    if args.positions_file is not None:
        positions_counter = PositionsCounter(args.positions_file, counter)    
    
    if args.chromosome is None:
        if args.positions_file is None:
            refs = counter.refs
        else:
            refs = positions_counter.refs
    else:
        refs = [args.chromosome, ]
        
    for ref in sorted(refs):
        positions_iter = counter.get_ref_iterator(ref)
        
        if args.positions_file is not None:
            positions_iter = positions_counter.get_ref_iterator(positions_iter)
        
        for row in positions_iter:
            data = row._data
            
            d_N = data._a_N + data._b_N
            d_T = data._a_T + data._b_T
            
            if  (min_d_N <= d_N <= max_d_N) and (min_d_T <= d_T <= max_d_T):
                i += 1
                
                if i % skip_size == 0:
                    data_set.append(data)
    
    return data_set

#=======================================================================================================================
# Helper factory classes for setting up parsers and classifiers.
#=======================================================================================================================
class ModelFactory(object):
    def get_counter(self, args):
        min_base_qual = args.min_base_qual
        min_map_qual = args.min_map_qual
        
        if args.model in ['binomial', 'beta_binomial']:
            qualities = 0            
        elif args.model in ['snvmix2', ]:
            qualities = 1
        
        normal_bam = BamFile(args.normal_file)
        tumour_bam = BamFile(args.tumour_file)
        ref_genome = FastaFile(args.reference_genome_file)
        
        
        counter = JointBinaryCounter(normal_bam,
                                     tumour_bam,
                                     ref_genome,
                                     min_base_qual,
                                     min_map_qual,
                                     qualities)
        
        return counter

    def get_model(self, args):
        priors = self._get_joint_snv_mix_priors(args)
        params = self._get_joint_snv_mix_params(args)
        
        if args.model == 'binomial':
            return BinomialModel(priors, params)
        elif args.model == 'beta_binomial':
            return BetaBinomialModel(priors, params)
        elif args.model == 'snvmix2':
            return SnvMixTwoModel(priors, params)
    
    def _get_joint_snv_mix_priors(self, args):
        if args.model in ['binomial', 'snvmix2']:
            priors = BinomialPriors()
        elif args.model in ['beta_binomial', ]:
            priors = BetaBinomialPriors()
        
        if args.mode == 'train' and args.priors_file is not None:
            priors.read_from_file(args.priors_file)
            
            sys.stderr.write("Using custom initial priors from {0}\n".format(args.priors_file))
        
        return priors
    
    def _get_joint_snv_mix_params(self, args):
        if args.model in ['binomial', 'snvmix2']:
            params = BinomialParameters()
        elif args.model in ['beta_binomial', ]:
            params = BetaBinomialParameters()    
        
        if args.mode == 'train' and args.initial_parameters_file is not None:            
            params.read_from_file(args.initial_parameters_file)
            
            sys.stderr.write("Using custom initial parameters from {0}\n".format(args.initial_parameters_file))
        
        if args.mode == 'classify' and args.parameters_file is not None:
            params.read_from_file(args.parameters_file)
            
            sys.stderr.write("Using custom parameters from {0}\n".format(args.parameters_file))
        
        return params
