'''
Created on 2011-03-31

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2

from joint_snv_mix.classification.base import *

from joint_snv_mix.classification.utils.data import JointData

class JointLatentVariables(EMLatentVariables):
    def _init_responsibilities(self, data):
        '''
        Intialise responsibilities via k-means clustering.
        '''
        a_1 = np.asarray(data.a['normal'], dtype=np.float64)
        b_1 = np.asarray(data.b['normal'], dtype=np.float64)
        p_1 = a_1 / (a_1 + b_1)
              
        a_2 = np.asarray(data.a['tumour'], dtype=np.float64)
        b_2 = np.asarray(data.b['tumour'], dtype=np.float64)
        p_2 = a_2 / (a_2 + b_2)

        shape = (data.nrows, 9)
        
        responsibilities = np.zeros(shape)
        
        init_centers = np.array((1., 0.5, 0.))
        
        cluster_centers_1, labels_1 = kmeans2(p_1, init_centers, minit='matrix')
        cluster_centers_2, labels_2 = kmeans2(p_2, init_centers, minit='matrix')

        labels = 3 * labels_1 + labels_2

        for id in range(9):
            index = labels == id
            
            responsibilities[index, id] = 1.
        
        self.responsibilities = responsibilities
        
class JointModelPriorParser(PriorParser):
    def __init__(self):
        PriorParser.__init__(self)
        
        self.ncomponent = self.nclass['normal'] * self.nclass['tumour']
    
    def _load_mix_weight_priors(self):       
        self.priors['kappa'] = np.zeros((self.ncomponent,))
            
        for i, genotype_tuple in enumerate(constants.joint_genotypes):
            genotype = "_".join(genotype_tuple)
        
            self.priors['kappa'][i] = self.parser.getfloat('kappa', genotype)
            
class JointParameterParser(ParameterParser):
    def __init__(self):
        ParameterParser.__init__(self)
        
        self.ncomponent = self.nclass['normal'] * self.nclass['tumour']

    def _load_mix_weights(self):       
        pi = np.zeros((self.ncomponent,))
            
        for i, genotype_tuple in enumerate(constants.joint_genotypes):
            genotype = "_".join(genotype_tuple)
        
            pi[i] = self.parser.getfloat('pi', genotype)
            
        self.parameters['pi'] = pi / pi.sum()            