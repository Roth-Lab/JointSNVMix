#!/usr/bin/env python

#=======================================================================================================================
# Joint SNVMix
# Author : Andrew Roth
#=======================================================================================================================

import argparse

from joint_snv_mix.classification.model_runners import run_classifier

from joint_snv_mix.pre_processing.mpileup_to_jcnt import main as mpileup_to_jcnt

from joint_snv_mix.pre_processing.mpileup_to_mcnt import main as mpileup_to_mcnt

from joint_snv_mix.post_processing.call_joint_genotypes import main as call_genotypes

from joint_snv_mix.post_processing.extract_jsm_positions import main as extract_jsm_positions

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

parser_mcnt.set_defaults( func=mpileup_to_mcnt )

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

parser_classify.add_argument( '--density', choices=['binomial', 'beta_binomial'], default='beta_binomial',
                              help='Density to be used in model.' )

parser_classify.add_argument( '--inference_algorithm', choices=['em', 'vb'], default='em',
                              help='Method used to infer model parameters.' )

train_group.set_defaults( func=run_classifier )

#===============================================================================
# Add call sub-command
#===============================================================================
parser_call = subparsers.add_parser( 'call',
                                     help="Call germlines, somatics or LOH sites based on user criteria." )

parser_call.add_argument( 'jsm_file_name',
                          help='Input JSM file name.' )

parser_call.add_argument( 'call_file_name',
                          help='Output file name.' )

parser_call.add_argument( '--genotype_class', choices=['Somatic', 'Germline', 'LOH'], default='Somatic',
                          help='Joint genotype call to call from file.' )

method_group = parser_call.add_mutually_exclusive_group( required=True )
method_group.add_argument( '-p', dest='prob_threshold', default=0.95,
                           help='Minimum probability to call a site in the given class.' )

method_group.add_argument( '-a', dest='argmax', default=False, action='store_true',
                           help='If set call call a site to a given class by argmax rule.' )

parser_call.set_defaults( func=call_genotypes )

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

if args.priors_file is None:
    args.train = False
else:
    args.train = True

args.func( args )
