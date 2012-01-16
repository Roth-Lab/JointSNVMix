#!/usr/bin/env python

#=======================================================================================================================
# JointSNVMix
# Author : Andrew Roth
#=======================================================================================================================
import argparse

from joint_snv_mix.runners.train import joint_snv_mix_one_train, joint_snv_mix_two_train
from joint_snv_mix.runners.classify import joint_snv_mix_one_classify, joint_snv_mix_two_classify 

parser = argparse.ArgumentParser(prog='JointSNVMix')
subparsers = parser.add_subparsers()

#=======================================================================================================================
# Parent parser for train and classifiy
#=======================================================================================================================
train_classify_parent = argparse.ArgumentParser(add_help=False)

train_classify_parent.add_argument('reference_genome',
                                   help='Path to reference genome fasta file.')

train_classify_parent.add_argument('normal_bam',
                                   help='''Path to BAM file for normal sample.''')

train_classify_parent.add_argument('tumour_bam',
                                   help='''Path to BAM file for tumour sample.''')

train_classify_parent.add_argument('--positions_file', default=None,
                                   help='''Path to a file containing a list of positions to create use for analysis.
                                   Should be space separated chrom pos. Additionally for each chromosome the positions
                                   should be sorted. The same format as samtools.''')

#=======================================================================================================================
# Parent parser for train
#=======================================================================================================================
train_parent = argparse.ArgumentParser(add_help=False, parents=[train_classify_parent])

train_parent.add_argument('parameter_file',
                          help='Path to file where the estimated parameters will be written to.')

train_parent.add_argument('--priors_file',
                          help='Path to a file with priors for the model parameters.')

train_parent.add_argument('--initial_parameter_file',
                          help='Path to a file with initial parameter values for the model.')

train_parent.add_argument('--min_normal_depth', default=10, type=int,
                          help='''Minimum depth of coverage in normal sample for a site to be eligible for use in
                          training set. Default 10''')

train_parent.add_argument('--min_tumour_depth', default=10, type=int,
                          help='''Minimum depth of coverage in tumour sample for a site to be eligible for use in
                          training set. Default 10''')

train_parent.add_argument('--max_normal_depth', default=100, type=int,
                          help='''Maximum depth of coverage in normal sample for a site to be eligible for use in
                          training set. Default 100''')

train_parent.add_argument('--max_tumour_depth', default=100, type=int,
                          help='''Maximum depth of coverage in tumour sample for a site to be eligible for use in
                          training set. Default 100''')

train_parent.add_argument('--max_iters', default=1000, type=int,
                          help='''Maximum number of iterations to used for training model. Default 1000''')

train_parent.add_argument('--skip_size', default=100, type=int,
                          help='''When subsampling will skip over this number of position before adding a site to the
                          subsample. Larger values lead to smaller subsample data sets with faster training and less
                          memory. Smaller values should lead to better parameter estimates. Default 100.''')

train_parent.add_argument('--convergence_threshold', default=1e-6, type=float,
                          help='''Convergence threshold for EM training. Once the change in objective function is below
                          this value training will end. Default 1e-6''')

#=======================================================================================================================
# Parent parser for classify command.
#=======================================================================================================================
classify_parent = argparse.ArgumentParser(add_help=False, parents=[train_classify_parent])

classify_parent.add_argument('out_file_name',
                             help='''Path to output file to be created. Output file will be tab separated file with
                             position information and class labels.''')

classify_parent.add_argument('parameter_file_name',
                              help='Path to file with the parameters to be used for classification.')

classify_parent.add_argument('out_file_name',
                              help='''Path to output file to be created. Output file will be tab separated file with
                              position information and class labels.''')

#=======================================================================================================================
# Parent parser for any classifier which uses simple count data
#=======================================================================================================================
count_parent = argparse.ArgumentParser(add_help=False)

count_parent.add_argument('--min_base_qual', default=10, type=int,
                          help='''Remove bases with base quality lower than this. Default is 10.''')

count_parent.add_argument('--min_map_qual', default=10, type=int,
                          help='''Remove bases with mapping quality lower than this. Default is 10.''')


#=======================================================================================================================
# Parsers for training subcommand
#=======================================================================================================================
train_parser = subparsers.add_parser('train', add_help=False,
                                    help='''Learn the parameters for a model from the SnvMix family. The output of this
                                    command can be passed as an option to some models under the classify command.''')

train_subparsers = train_parser.add_subparsers()

#---------------------------------------------------------------------------------------------------------------------- 
train_jsm1_parser = train_subparsers.add_parser('joint_snv_mix_one', parents=[train_parent, count_parent],
                                                help='''Learn the model parameters for the JointSnvMix1 model.''')

train_jsm1_parser.set_defaults(func=joint_snv_mix_one_train)

#---------------------------------------------------------------------------------------------------------------------- 
train_jsm2_parser = train_subparsers.add_parser('joint_snv_mix_two', parents=[train_parent],
                                                help='''Learn the model parameters for the JointSnvMix2 model.''')

train_jsm2_parser.set_defaults(func=joint_snv_mix_two_train)

#=======================================================================================================================
# Parsers for classify subcommand
#=======================================================================================================================
classify_parser = subparsers.add_parser('classify',
                                        help='''Classify the data from two BAM files from normal/tumour pair using one
                                        of the several strategies.''')

classify_subparsers = classify_parser.add_subparsers()

#---------------------------------------------------------------------------------------------------------------------- 
classify_jsm1_parser = classify_subparsers.add_parser('joint_snv_mix_one', parents=[classify_parent, count_parent],
                                                      help='''Classify paired data using the JointSnvMix1 model.''')

classify_jsm1_parser.set_defaults(func=joint_snv_mix_one_classify)

#---------------------------------------------------------------------------------------------------------------------- 
classify_jsm2_parser = classify_subparsers.add_parser('joint_snv_mix_two', parents=[classify_parent],
                                                      help='''Classify paired data using the JointSnvMix2 model.''')

classify_jsm2_parser.set_defaults(func=joint_snv_mix_two_classify)

#=======================================================================================================================
# Parse arguments.
#=======================================================================================================================
args = parser.parse_args()

args.func(args)
