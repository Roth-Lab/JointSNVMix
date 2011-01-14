'''
Created on 2011-01-13

@author: Andrew Roth
'''
from bin.independent_beta_binomial_em import main as run_independent_beta_binomial

from bin.independent_binomial_em import main as run_independent_binomial

from bin.joint_binomial_em import main as run_joint_binomial

from bin.joint_beta_binomial_em import main as run_joint_beta_binomial

from bin.joint_binomial_vb import main as run_joint_binomial_vb

from bin.beta_binomial_mixture_em import main as run_beta_binomial_mixture

def run_classifier( args ):
    if args.model == "independent":
        if args.density == "binomial":
            run_independent_binomial( args )
        elif args.density == "beta_binomial":
            run_independent_beta_binomial( args )
            
    elif args.model == "joint":
        if args.density == "binomial":
            if args.inference_algorithm == "vb":
                run_joint_binomial_vb( args )
            else:                
                run_joint_binomial( args )
                
        elif args.density == "beta_binomial":
            run_joint_beta_binomial( args )
    
    elif args.model == "mixture":
        if args.density == "beta_binomial":
            run_beta_binomial_mixture( args )
