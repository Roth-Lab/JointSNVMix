'''
Created on 2011-01-31

@author: andrew
'''
import numpy as np
from joint_snv_mix.classification.random_variable import SixParameterBetaBinomialMixture
from joint_snv_mix.classification.data import JointData
from tests.snp_simulators.beta_binomial import draw_bb_sample, draw_joint_sample

priors_file = "/home/andrew/workspace/joint_snv_mix/config/joint_bb_em.priors.cfg"

mm = SixParameterBetaBinomialMixture( priors_file=priors_file )

params = {}
params['alpha'] = np.array( [
                             [99, 5, 1],
                             [99, 5, 1]
                             ], dtype=np.float )

params['beta'] = np.array( [
                             [1, 5, 99],
                             [1, 5, 99]
                             ], dtype=np.float )

params['pi'] = np.array( [1e5, 100, 10, 100, 1e3, 100, 1, 5, 1e3], dtype=np.float )
params['pi'] = params['pi'] / params['pi'].sum()

data, labels = draw_joint_sample( params, 100000, 40 )

data = JointData( data )

mm.train( data )

pred_labels = np.argmax( mm.classify( data ) )

print pred_labels == labels
mm.print_parameters()