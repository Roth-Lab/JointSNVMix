'''
Created on 2011-08-11

@author: Andrew Roth
'''
import pysam

from joint_snv_mix.counters.joint_binary_counter import JointBinaryBaseCounter
from joint_snv_mix.counters.joint_binary_quality_counter import JointBinaryQualityCounter
from joint_snv_mix.counters.positions_counter import PositionsCounter

from joint_snv_mix.trainers.joint_snv_mix import JointSnvMixOneModel, JointSnvMixTwoModel, JointSnvMixOneSubsampler, \
    JointSnvMixTwoSubsampler, JointSnvMixParameters, JointSnvMixPriors

from joint_snv_mix.trainers.snv_mix import PairedSnvMixOneModel, PairedSnvMixTwoModel, SnvMixOneSubsampler, \
    SnvMixTwoSubsampler, PairedSnvMixParameters, PairedSnvMixPriors

def snv_mix_one_train(args):
    counter = get_base_counter(args)
    ss = SnvMixOneSubsampler(args.skip_size, args.min_normal_depth, args.min_tumour_depth)
    sample = ss.subsample(counter)
    
    params = get_indep_params(args)
    
    model = PairedSnvMixOneModel(params)
    
    train(model, sample, args)
    
def snv_mix_two_train(args): 
    counter = get_qual_counter(args)
    ss = SnvMixTwoSubsampler(args.skip_size, args.min_normal_depth, args.min_tumour_depth)
    sample = ss.subsample(counter)
    
    params = get_indep_params(args)
    
    model = PairedSnvMixTwoModel(params)
    
    train(model, sample, args)

def joint_snv_mix_one_train(args): 
    counter = get_base_counter(args)
    ss = JointSnvMixOneSubsampler(args.skip_size, args.min_normal_depth, args.min_tumour_depth)
    sample = ss.subsample(counter)
    
    params = get_joint_params(args)
    model = JointSnvMixOneModel(params)
    
    train(model, sample, args)

    
def joint_snv_mix_two_train(args): 
    counter = get_qual_counter(args)
    
    params = get_joint_params(args)
    
    model = JointSnvMixTwoModel(params)
    
    ss = JointSnvMixTwoSubsampler(args.skip_size, args.min_normal_depth, args.min_tumour_depth)
    sample = ss.subsample(counter)
    
    train(model, sample, args)

def get_base_counter(args):
    normal_bam = pysam.Samfile(args.normal_bam, 'rb')
    tumour_bam = pysam.Samfile(args.tumour_bam, 'rb')
    ref_genome = pysam.Fastafile(args.reference_genome)
    
    counter = JointBinaryBaseCounter(normal_bam,
                                     tumour_bam,
                                     ref_genome,
                                     min_base_qual=args.min_base_qual,
                                     min_map_qual=args.min_map_qual)

    if args.positions_file is not None:
        counter = PositionsCounter(counter)
    
    return counter
    
def get_qual_counter(args):
    normal_bam = pysam.Samfile(args.normal_bam, 'rb')
    tumour_bam = pysam.Samfile(args.tumour_bam, 'rb')
    ref_genome = pysam.Fastafile(args.reference_genome)
    
    counter = JointBinaryQualityCounter(normal_bam, tumour_bam, ref_genome)
    
    if args.positions_file is not None:
        counter = PositionsCounter(counter)
    
    return counter

def get_indep_params(args):
    priors = PairedSnvMixPriors()
    priors.read_from_file(args.priors_file_name)
    
    params = PairedSnvMixParameters(priors=priors)
    params.read_from_file(args.initial_parameter_file_name)
    
    return params

def get_joint_params(args):
    priors = JointSnvMixPriors()
    priors.read_from_file(args.priors_file_name)
    
    params = JointSnvMixParameters(priors=priors)
    params.read_from_file(args.initial_parameter_file_name)
    
    return params

def train(model, sample, args):        
    model.fit(sample, args.convergence_threshold, args.max_iters)

    model.params.write_to_file(args.parameter_file_name)
    
