import numpy as np

from joint_snv_mix.classification.classification import JointBetaBinomialModel, parse_joint_priors_file
from tests.simulation.beta_binomial import draw_joint_sample
from joint_snv_mix.classification.data import JointData

priors_file = "/home/andrew/workspace/joint_snv_mix/config/joint_bb_em.priors.cfg"

density = "beta_binomial"
n = int( 1e5 )
mean_depth = 40

alpha = {}
alpha['normal'] = np.array( [99, 1000, 1 ], dtype=np.float )
alpha['tumour'] = np.array( [99, 30, 1 ], dtype=np.float )

beta = {}
beta['normal'] = np.array( [0.3, 90, 99], dtype=np.float )
beta['tumour'] = np.array( [1, 20, 10], dtype=np.float )

pi = np.array( [1e4, 100, 100,
                100, 1e3, 100,
                1, 10, 1e3] )

pi = pi / pi.sum()

priors = parse_joint_priors_file( priors_file, n, density )

model = JointBetaBinomialModel()

counts, labels = draw_joint_sample( alpha, beta, pi, n, mean_depth )
data = JointData( counts )

params = model.train( data, priors, 1000, 1e-6 )

print params
