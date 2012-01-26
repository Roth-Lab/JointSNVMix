#!/usr/bin/env python

#=======================================================================================================================
# JointSNVMix
# Author : Andrew Roth
#=======================================================================================================================
import argparse

from joint_snv_mix.run import classify, train 

parser = argparse.ArgumentParser(prog='JointSNVMix')
subparsers = parser.add_subparsers()

#=======================================================================================================================
# Parent parser for train and classifiy
#=======================================================================================================================
train_classify_parent = argparse.ArgumentParser(add_help=False)

train_classify_parent.add_argument('reference_genome_file',
                                   help='Path to reference genome fasta file.')

train_classify_parent.add_argument('normal_file',
                                   help='''Path to BAM file for normal sample.''')

train_classify_parent.add_argument('tumour_file',
                                   help='''Path to BAM file for tumour sample.''')

train_classify_parent.add_argument('--chromosome', default=None,
                                   help='''Chromosome to analyse. If not set all chromosomes will be analysed.''')

train_classify_parent.add_argument('--min_base_qual', default=0, type=int,
                                   help='''Remove bases with base quality lower than this. Default is 0.''')

train_classify_parent.add_argument('--min_map_qual', default=0, type=int,
                                   help='''Remove bases with mapping quality lower than this. Default is 0.''')

train_classify_parent.add_argument('--model', default='snvmix2', choices=['binomial', 'beta_binomial', 'snvmix2'],
                                   help='''Density to be used in mixture model. Default is snvmix2.''')

train_classify_parent.add_argument('--positions_file', default=None,
                                   help='''Path to a file containing a list of positions to create use for analysis.
                                   Should be space separated chrom pos. Additionally for each chromosome the positions
                                   should be sorted. The same format as samtools.''')

#=======================================================================================================================
# Parent parser for train
#=======================================================================================================================
train_parser = subparsers.add_parser('train', parents=[train_classify_parent],
                                    help='''Learn the parameters for a model from the SnvMix family. The output of this
                                    command can be passed as an option to some models under the classify command.''')

train_parser.add_argument('estimated_parameters_file',
                          help='Path to file where the estimated parameters will be written to.')

train_parser.add_argument('--priors_file',
                          help='Path to a file with priors for the model parameters.')

train_parser.add_argument('--initial_parameters_file',
                          help='Path to a file with initial parameter values for the model.')

train_parser.add_argument('--min_normal_depth', default=10, type=int,
                          help='''Minimum depth of coverage in normal sample for a site to be eligible for use in
                          training set. Default 10''')

train_parser.add_argument('--min_tumour_depth', default=10, type=int,
                          help='''Minimum depth of coverage in tumour sample for a site to be eligible for use in
                          training set. Default 10''')

train_parser.add_argument('--max_normal_depth', default=100, type=int,
                          help='''Maximum depth of coverage in normal sample for a site to be eligible for use in
                          training set. Default 100''')

train_parser.add_argument('--max_tumour_depth', default=100, type=int,
                          help='''Maximum depth of coverage in tumour sample for a site to be eligible for use in
                          training set. Default 100''')

train_parser.add_argument('--max_iters', default=1000, type=int,
                          help='''Maximum number of iterations to used for training model. Default 1000''')

train_parser.add_argument('--skip_size', default=1, type=int,
                          help='''When subsampling will skip over this number of position before adding a site to the
                          subsample. Larger values lead to smaller subsample data sets with faster training and less
                          memory. Smaller values should lead to better parameter estimates. Default 1.''')

train_parser.add_argument('--convergence_threshold', default=1e-6, type=float,
                          help='''Convergence threshold for EM training. Once the change in objective function is below
                          this value training will end. Default 1e-6''')

train_parser.set_defaults(func=train, mode='train')

#=======================================================================================================================
# Parent parser for classify command.
#=======================================================================================================================
classify_parser = subparsers.add_parser('classify', parents=[train_classify_parent],
                                        help='''Classify the data from two BAM files from normal/tumour pair using one
                                        of the several strategies.''')

classify_parser.add_argument('--out_file',
                             help='''If set results will be output here, otherwise they will go to STDOUT. Output file
                             will be tab separated file with position information and class labels.''')

classify_parser.add_argument('--print_all_positions', action='store_true', default=False,
                             help='''By default only positions with a variant read in the tumour are printed. If this
                             flag is set all sites will be printed.''')

classify_parser.add_argument('--parameters_file',
                             help='Path to a file with custom parameters values for the model.')

classify_parser.add_argument('--somatic_threshold', default=0.0, type=float,
                             help='''Only sites with P(Somatic) = p_AA_AB + p_AA_BB greater than equal this value
                             will be printed. Default is 0.''' )

classify_parser.set_defaults(func=classify, mode='classify')

#=======================================================================================================================
# Parse arguments.
#=======================================================================================================================
args = parser.parse_args()

args.func(args)
