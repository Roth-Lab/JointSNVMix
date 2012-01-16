'''
Created on 2011-08-14

@author: Andrew Roth
'''
import pysam

from joint_snv_mix.counters.joint_binary_counter import JointBinaryBaseCounter
from joint_snv_mix.counters.joint_binary_quality_counter import JointBinaryQualityCounter
from joint_snv_mix.counters.positions_counter import PositionsCounter

from joint_snv_mix.classifiers.classifier import classify_counter
from joint_snv_mix.classifiers.joint_snv_mix import JointSnvMixOneClassifier
from joint_snv_mix.classifiers.joint_snv_mix_qualities import JointSnvMixTwoClassifier

from joint_snv_mix.trainers.snv_mix import PairedSnvMixParameters
from joint_snv_mix.trainers.joint_snv_mix import JointSnvMixParameters
   
def joint_snv_mix_one_classify(args):
    counter = get_base_counter(args)
    
    params = JointSnvMixParameters()
    params.read_from_file(args.parameter_file_name)
    
    classifier = JointSnvMixOneClassifier(mu_N=params.mu_N,
                                          mu_T=params.mu_T,
                                          pi=params.pi)
    
    classify(counter, classifier, args)
    
def joint_snv_mix_two_classify(args):
    counter = get_qual_counter(args)
    
    params = JointSnvMixParameters()
    params.read_from_file(args.parameter_file_name)
    
    classifier = JointSnvMixTwoClassifier(mu_N=params.mu_N,
                                          mu_T=params.mu_T,
                                          pi=params.pi)
    
    classify(counter, classifier, args)
    
#=======================================================================================================================
# Helper functions
#=======================================================================================================================
def classify(counter, classifier, args):
    classify_counter(counter, classifier, args.out_file_name)
    
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
        counter = PositionsCounter(counter, args.positions_file)
    
    return counter
    
def get_qual_counter(args):
    normal_bam = pysam.Samfile(args.normal_bam, 'rb')
    tumour_bam = pysam.Samfile(args.tumour_bam, 'rb')
    ref_genome = pysam.Fastafile(args.reference_genome)
    
    counter = JointBinaryQualityCounter(normal_bam, tumour_bam, ref_genome)
    
    if args.positions_file is not None:
        counter = PositionsCounter(counter, args.positions_file)
    
    return counter
