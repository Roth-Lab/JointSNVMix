'''
Created on 2011-03-31

@author: Andrew Roth
'''
import numpy as np

from scipy.cluster.vq import kmeans2

from joint_snv_mix import constants

from joint_snv_mix.classification.base import *

from joint_snv_mix.classification.utils.data import IndependentData

class PairedIndependentModel(object):
    def __init__(self, model):
        self.model = model
        
    def train(self, joint_data, priors, max_iters, tolerance):
        parameters = {}
                
        for genome in constants.genomes:
            data = joint_data.get_independent_data(genome)
            
            parameters[genome] = self.model.train(data, priors[genome], max_iters, tolerance)
        
        return parameters
            
    def classify(self, joint_data, parameters):
        self.resp = {}
        
        for genome in constants.genomes:
            data = joint_data.get_independent_data(genome)
            
            self.resp[genome] = self.model.classify(data, parameters[genome])
        
        return self._get_joint_resp()
        
    def _get_joint_resp(self):
        joint_resp = []
        
        nrows = self.resp['normal'].shape[0]
        nclass = self.resp['normal'].shape[1]
        
        col_shape = (nrows, 1)
        
        for i in range(nclass):
            joint_resp.append(self.resp['normal'][:, i].reshape(col_shape) + self.resp['tumour'])
        
        joint_resp = np.hstack(joint_resp)
        
        return joint_resp
            
class IndependentLatenVariables(EMLatentVariables):        
    def _init_responsibilities(self, data):
        a = np.asarray(data.a, dtype=np.float64)
        b = np.asarray(data.b, dtype=np.float64)
        
        shape = (data.nrows, 3)
        
        responsibilities = np.zeros(shape)
        
        p = a / (a + b)
        
        init_centers = np.array([1., 0.5, 0.])
        
        cluster_centers, labels = kmeans2(p, init_centers)
        
        sorted_centers = np.argsort(cluster_centers)
        
        for id in sorted_centers:
            index = labels == id
            
            responsibilities[index, id] = 1.0
        
        self.responsibilities = responsibilities

class IndepedendentPriorParser(PriorParser):    
    def _load_mix_weight_priors(self):       
        for genome in constants.genomes:            
            self.priors[genome]['kappa'] = np.zeros((self.nclass[genome],))
            
            for i, genotype in enumerate(constants.genotypes):                
                genome_genotype = "_".join((genome, genotype))
            
                self.priors[genome]['kappa'][i] = self.parser.getfloat('kappa', genome_genotype)
                
class IndepedendentParameterParser(ParameterParser):    
    def _load_mix_weights(self):       
        for genome in constants.genomes:            
            pi = np.zeros((self.nclass[genome],))
            
            for i, genotype in enumerate(constants.genotypes):                
                genome_genotype = "_".join((genome, genotype))
            
                pi[i] = self.parser.getfloat('pi', genome_genotype)
                
            self.parameters[genome]['pi'] = pi / pi.sum()
