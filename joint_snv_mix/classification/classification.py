'''
Created on 2011-01-13

@author: Andrew Roth
'''
from bin.independent_beta_binomial_em import main as run_independent_beta_binomial

from bin.independent_binomial_em import main as run_independent_binomial

from bin.joint_binomial_em import main as run_joint_binomial

from bin.joint_beta_binomial_em import main as run_joint_beta_binomial

def run_classifier( args ):
    if args.model == "independent":
        if args.density == "binomial":
            run_independent_binomial( args )
        elif args.density == "beta_binomial":
            run_independent_beta_binomial( args )
            
    elif args.model == "joint":
        if args.density == "binomial":
            run_joint_binomial( args )
        elif args.density == "beta_binomial":
            run_joint_beta_binomial( args )
