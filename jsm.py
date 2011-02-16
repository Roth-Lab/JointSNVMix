#!/usr/bin/env python

#=======================================================================================================================
# Joint SNVMix
# Author : Andrew Roth
#=======================================================================================================================

import argparse

from joint_snv_mix.classification.model_runners import run_classifier

from joint_snv_mix.pre_processing.mpileup_to_jcnt import main as mpileup_to_jcnt

from joint_snv_mix.pre_processing.mpileup_to_mcnt import main as mpileup_to_mcnt

from joint_snv_mix.pre_processing.varscan_to_jcnt import main as varscan_to_jcnt

from joint_snv_mix.post_processing.call_somatics import main as call_somatics

from joint_snv_mix.post_processing.extract_jsm_positions import main as extract_jsm_positions
from joint_snv_mix.classification.conan import run_conan

parser = argparse.ArgumentParser( prog='JointSNVMix' )
subparsers = parser.add_subparsers()

#===============================================================================
# Add jcnt sub-command
#===============================================================================
parser_jcnt = subparsers.add_parser( 'jcnt',
                                     help='Convert mpileup files to jcnt file format.' )

parser_jcnt.add_argument( 'mpileup_file_name',
                          help='''Samtools mpileup format file. When creating file with samtools the first bam file
                          passed as arguments should be the normal and the second the tumour file.''' )

parser_jcnt.add_argument( 'jcnt_file_name',
                          help='Name of joint counts (jcnt) output files to be created.' )

parser_jcnt.add_argument( '--min_depth', default=1, type=int,
                          help='''Minimum depth of coverage in both tumour and normal sample required to use a site in
                          the analysis.''' )

parser_jcnt.add_argument( '--min_qual', default=13, type=int,
                          help='''Remove bases with base qualities lower than this value. Note if samtools calmd is used
                          to pre-process the bam then BAQ qualities may replace the base qualities in the 
                          mpileup file.''' )

parser_jcnt.add_argument( '--bzip2', action='store_true',
                          help='''Set if file is in bzip2 format.''' )

parser_jcnt.set_defaults( func=mpileup_to_jcnt )

#===============================================================================
# Add mcnt sub-command
#===============================================================================
parser_mcnt = subparsers.add_parser( 'mcnt',
                                     help='Convert mpileup files to mcnt file format.' )

parser_mcnt.add_argument( 'mpileup_file_name',
                          help='''Samtools mpileup format file. When creating file with samtools the first bam file
                          passed as arguments should be the normal and the second the tumour file.''' )

parser_mcnt.add_argument( 'mcnt_file_name',
                          help='Name of joint counts (jcnt) output files to be created.' )

parser_mcnt.add_argument( '--min_depth', default=1, type=int,
                          help='''Minimum depth of coverage in both tumour and normal sample required to use a site in
                          the analysis.''' )

parser_mcnt.add_argument( '--bzip2', action='store_true',
                          help='''Set if file is in bzip2 format.''' )

parser_mcnt.set_defaults( func=mpileup_to_mcnt )

#===============================================================================
# Add varscan sub-command
#===============================================================================
parser_mcnt = subparsers.add_parser( 'varscan',
                                     help='Convert varscan files to jcnt file format.' )

parser_mcnt.add_argument( 'varscan_file_name',
                          help='''Samtools mpileup format file. When creating file with samtools the first bam file
                          passed as arguments should be the normal and the second the tumour file.''' )

parser_mcnt.add_argument( 'jcnt_file_name',
                          help='Name of joint counts (jcnt) output files to be created.' )

parser_mcnt.set_defaults( func=varscan_to_jcnt )

#===============================================================================
# Add analyse sub-command
#===============================================================================
parser_classify = subparsers.add_parser( 'classify',
                                        help='''Run a JointSNVMix analysis. Requires that a jcnt file has been
                                        created''' )

parser_classify.add_argument( 'jcnt_file_name',
                             help='Name of joint counts (jcnt) file to be used as input.' )

parser_classify.add_argument( 'jsm_file_name',
                             help='Name of JointSNVMix (jsm) output files to be created.' )

file_group = parser_classify.add_mutually_exclusive_group( required=True )

file_group.add_argument( '--params_file', default=None,
                         help='''File containing model parameters to use for classification.
                         If set the model will not be trained.''' )

file_group.add_argument( '--priors_file', default=None,
                         help='File containing prior distribution parameters to use for training. \
                         If set the model will be trained.' )

train_group = parser_classify.add_argument_group( title='Training Parameters',
                                                 description='Options for training the model.' )

train_group.add_argument( '--max_iters', default=1000, type=int,
                          help='''Maximum number of iterations to used for training model. Default 1000''' )

train_group.add_argument( '--subsample_size', default=0, type=int,
                          help='''Size of random subsample to use for training. If not set the whole data set will be
                          used.''' )

train_group.add_argument( '--convergence_threshold', default=1e-6, type=float,
                          help='''Convergence threshold for EM training. Once the change in objective function is below
                          this value training will end. Defaul 1e-6''' )

parser_classify.add_argument( '--model', choices=['independent', 'joint', 'chromosome'], default='joint',
                              help='Model type to use for classification.' )

parser_classify.add_argument( '--density', choices=['binomial', 'beta_binomial', 'multinomial'], default='beta_binomial',
                              help='Density to be used in model.' )

parser_classify.add_argument( '--inference_algorithm', choices=['em', 'vb'], default='em',
                              help='Method used to infer model parameters.' )

train_group.set_defaults( func=run_classifier )

#===============================================================================
# Add conan sub-command
#===============================================================================
parser_conan = subparsers.add_parser( 'conan',
                                        help='''Run a ConanSNVMix analysis. Requires that a cncnt file has been
                                        created''' )

parser_conan.add_argument( 'cncnt_file_name',
                             help='Name of conan counts (cncnt) file to be used as input.' )

parser_conan.add_argument( 'cnsm_file_name',
                             help='Name of ConanSNVMix (cnsm) output files to be created.' )

parser_conan.add_argument( '--density', choices=['binomial', 'beta_binomial'], default='beta_binomial',
                              help='Density to be used in model.' )

train_group = parser_conan.add_argument_group( title='Training Parameters',
                                                 description='Options for training the model.' )

train_group.add_argument( '--max_iters', default=1000, type=int,
                          help='''Maximum number of iterations to used for training model. Default 1000''' )

train_group.add_argument( '--subsample_size', default=0, type=int,
                          help='''Size of random subsample to use for training. If not set the whole data set will be
                          used.''' )

train_group.add_argument( '--convergence_threshold', default=1e-6, type=float,
                          help='''Convergence threshold for EM training. Once the change in objective function is below
                          this value training will end. Defaul 1e-6''' )

train_group.set_defaults( func=run_conan )

#===============================================================================
# Add call sub-command
#===============================================================================
parser_call = subparsers.add_parser( 'call',
                                     help="Call germlines, somatics or LOH sites based on user criteria." )

parser_call.add_argument( 'jsm_file_name',
                          help='Input JSM file name.' )

parser_call.add_argument( 'out_file_name',
                          help='Output file name.' )

parser_call.set_defaults( func=call_somatics )

#===============================================================================
# Add extract sub-command
#===============================================================================
parser_extract = subparsers.add_parser( 'extract',
                                        help='Extract a set of positions from a jsm file.' )

parser_extract.add_argument( 'jsm_file_name',
                             help='JSM file to extract positions from.' )

parser_extract.add_argument( 'positions_file_name',
                             help='''List of positions to extract. Format is tab delimited with the chromosome in the
                             first column and second position in the second i.e. "X    12345" ''' )

parser_extract.set_defaults( func=extract_jsm_positions )

#===============================================================================
# Run
#===============================================================================
args = parser.parse_args()

args.func( args )
