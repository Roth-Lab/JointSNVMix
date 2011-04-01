#!/usr/bin/env python

#=======================================================================================================================
# Joint SNVMix
# Author : Andrew Roth
#=======================================================================================================================
import argparse

from joint_snv_mix.classification.model_runners import run_snvmix

from joint_snv_mix.pre_processing.bam_to_jcnt import bam_to_jcnt

from joint_snv_mix.classification.threhsold_classifier import run_threshold

from joint_snv_mix.classification.fisher_classifier import run_fisher

parser = argparse.ArgumentParser(prog='JointSNVMix')
subparsers = parser.add_subparsers()

#===============================================================================
# Add jcnt sub-command
#===============================================================================
parser_jcnt = subparsers.add_parser('jcnt',
                                    help='Convert bam files to jcnt file format.')

parser_jcnt.add_argument('reference_genome_file_name',
                          help='Path to reference genome fasta.')

parser_jcnt.add_argument('normal_bam_file_name',
                          help='''Normal BAM file.''')

parser_jcnt.add_argument('tumour_bam_file_name',
                          help='''Tumour BAM file.''')

parser_jcnt.add_argument('jcnt_file_name',
                          help='Name of joint counts (jcnt) output files to be created.')

parser_jcnt.add_argument('--min_depth', default=1, type=int,
                          help='''Minimum depth of coverage in both tumour and normal sample required to use a site in
                          the analysis. Default is 1.''')

parser_jcnt.add_argument('--min_base_qual', default=10, type=int,
                          help='''Remove bases with base quality lower than this. Default is 10.''')

parser_jcnt.add_argument('--min_map_qual', default=10, type=int,
                          help='''Remove bases with mapping quality lower than this. Default is 10.''')

parser_jcnt.add_argument('--positions_file', default=None,
                          help='''Path to list of positions to create jcnt on. Should be space separated chrom pos. Same
                          as samtools.''')

parser_jcnt.set_defaults(func=bam_to_jcnt)


#===============================================================================
# Add snvmix model sub-command
#===============================================================================
parser_snvmix = subparsers.add_parser('snvmix',
                                        help='''Run a JointSNVMix or paired independent SNVMix analysis. Requires that a
                                        jcnt file has been created.''')

parser_snvmix.add_argument('jcnt_file_name',
                             help='Name of joint counts (jcnt) file to be used as input.')

parser_snvmix.add_argument('jsm_file_name',
                             help='Name of JointSNVMix (jsm) output files to be created.')

file_group = parser_snvmix.add_mutually_exclusive_group(required=True)

file_group.add_argument('--params_file', default=None,
                         help='''File containing model parameters to use for classification.
                         If set the model will not be trained.''')

file_group.add_argument('--priors_file', default=None,
                         help='File containing prior distribution parameters to use for training. \
                         If set the model will be trained.')

train_group = parser_snvmix.add_argument_group(title='Training Parameters',
                                                 description='Options for training the model.')

train_group.add_argument('--max_iters', default=1000, type=int,
                          help='''Maximum number of iterations to used for training model. Default 1000''')

train_group.add_argument('--subsample_size', default=0, type=int,
                          help='''Size of random subsample to use for training. If not set the whole data set will be
                          used.''')

train_group.add_argument('--convergence_threshold', default=1e-6, type=float,
                          help='''Convergence threshold for EM training. Once the change in objective function is below
                          this value training will end. Defaul 1e-6''')

parser_snvmix.add_argument('--model', choices=['independent', 'joint', 'chromosome'],
                              default='joint', help='Model type to use for classification.')

parser_snvmix.add_argument('--density', choices=['binomial', 'beta_binomial'], default='beta_binomial',
                              help='Density to be used in model.')

parser_snvmix.set_defaults(func=run_snvmix)

#===============================================================================
# Add fisher model sub-command
#===============================================================================
parser_fisher = subparsers.add_parser('fisher',
                                        help='''Run a fisher exact test based analysis. Requires that a jcnt file has
                                        been created.''')

parser_fisher.add_argument('jcnt_file_name',
                             help='Name of joint counts (jcnt) file to be used as input.')

parser_fisher.add_argument('tsv_file_name',
                             help='Name of output tsv file to write predictions to.')

parser_fisher.add_argument('--model', choices=['independent', 'joint'],
                              default='joint', help='Model type to use for classification.')

parser_fisher.add_argument('--p_value_threshold', default=0.05, type=float,
                              help='''Threshold for declaring a site homozygous in tumour.''')

parser_fisher.add_argument('--base_line_error', default=0.001, type=float,
                              help='''Expected error rate.''')

parser_fisher.add_argument('--min_var_freq', default=0.1, type=float,
                              help='''Minimum frequency for declaring a site a variant.''')

parser_fisher.add_argument('--min_hom_freq', default=0.8, type=float,
                              help='''Minimum frequency to decalare genotype BB.''')

parser_fisher.add_argument('--min_var_depth', default=4, type=int,
                              help='''Sites with fewer variant reads in the tumour than this will always be called
                              reference.''')

parser_fisher.set_defaults(func=run_fisher)

#===============================================================================
# Add threshold model sub-command
#===============================================================================
parser_threshold = subparsers.add_parser('threshold',
                                        help='''Run a naive threshold based classifier. Requires that a jcnt file has
                                        been created.''')

parser_threshold.add_argument('jcnt_file_name',
                              help='Name of joint counts (jcnt) file to be used as input.')

parser_threshold.add_argument('tsv_file_name',
                              help='Name of output tsv file to write predictions to.')

parser_threshold.add_argument('--normal_threshold', default=0.05, type=float,
                              help='''Threshold for declaring a site homozygous in normal.''')

parser_threshold.add_argument('--tumour_threshold', default=0.1, type=float,
                              help='''Threshold for declaring a site homozygous in tumour.''')

parser_threshold.add_argument('--min_var_depth', default=4, type=int,
                              help='''Sites with fewer variant reads in the tumour than this will always be called
                              reference.''')

parser_threshold.set_defaults(func=run_threshold)

#===============================================================================
# Add call_jsm sub-command
#===============================================================================
parser_call_jsm = subparsers.add_parser('call_somatics',
                                     help="Call somatics from a jsm file.")

parser_call_jsm.add_argument('jsm_file_name',
                          help='Input JSM file name.')

parser_call_jsm.add_argument('out_file_name',
                          help='Output file name.')

parser_call_jsm.add_argument('--threshold', default=0.99, type=float,
                          help='''Probability threshold for calling mutation somatics. Ignored if auto flag is set.''')

parser_call_jsm.add_argument('--auto', action='store_true', default=False,
                          help='''Automatically determine probability threshold to use for calling somatics.''')

parser_call_jsm.set_defaults(func=call_somatics_from_jsm)

#===============================================================================
# Add extract_positions sub-command
#===============================================================================
parser_extract_positions = subparsers.add_parser('extract_positions',
                                        help='Extract a set of positions from a jsm file and print to stdout.')

parser_extract_positions.add_argument('jsm_file_name',
                                      help='JSM file to extract positions from.')

parser_extract_positions.add_argument('positions_file_name',
                                      help='''List of positions to extract. Format is tab delimited with the chromosome
                                      in the first column and second position in the second i.e. "X    12345" ''')

parser_extract_positions.set_defaults(func=extract_jsm_positions)

#=======================================================================================================================
# Add extract_parameters sub_command
#=======================================================================================================================
parser_extract = subparsers.add_parser('extract_parameters',
                                        help='Extract a parameters from jsm file and print to stdout.')

parser_extract.add_argument('jsm_file_name',
                             help='JSM file to extract positions from.')

parser_extract.set_defaults(func=extract_jsm_parameters)

#===============================================================================
# Run
#===============================================================================
args = parser.parse_args()

args.func(args)
