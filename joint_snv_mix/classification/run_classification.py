'''
Created on 2011-04-04

@author: andrew
'''
from joint_snv_mix.classification.deterministic import ThresholdModel, IndependentFisherModel, DeterministicRunner, \
    JointFisherModel

from joint_snv_mix.classification.joint_binomial import JointBinomialModel, JointBinomialParameterParser, \
    JointBinomialPriorParser

from joint_snv_mix.classification.indep import PairedIndependentModel

from joint_snv_mix.classification.indep_binomial import IndependentBinomialModel, IndependentBinomialParameterParser, \
    IndependentBinomialPriorParser
from joint_snv_mix.classification.base import ProbabilisticModelRunner

def run_threshold(args):
    model = ThresholdModel(
                           args.normal_threshold,
                           args.tumour_threshold,
                           args.min_var_depth
                           )

    runner = DeterministicRunner(model)    
    runner.run(args.jcnt_file_name, args.tsv_file_name)

def run_fisher(args):   
    if args.model == "independent":
        model_class = IndependentFisherModel
    elif args.model == "joint":
        model_class = JointFisherModel
        
    model = model_class(
                        args.p_value_threshold,
                        args.base_line_error,
                        args.min_var_freq,
                        args.min_hom_freq,
                        args.min_var_depth
                        ) 
    
    runner = DeterministicRunner(model)
    runner.run(args.jcnt_file_name, args.tsv_file_name)

def run_binomial(args):
    if args.priors_file is not None:
        args.train = True
    else:
        args.train = False
        
    if args.model == "independent":
        model = PairedIndependentModel(IndependentBinomialModel())
        parameter_parser = IndependentBinomialParameterParser()
        priors_parser = IndependentBinomialPriorParser()
        
    elif args.model == "joint":
        model = JointBinomialModel()
        parameter_parser = JointBinomialParameterParser()
        priors_parser = JointBinomialPriorParser()
        
    runner = ProbabilisticModelRunner(model, priors_parser, parameter_parser)
    runner.run(args)
