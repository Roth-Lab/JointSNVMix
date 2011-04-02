'''
Created on 2010-11-19

@author: Andrew Roth
'''
class JointData(object):
    def __init__(self, X):
        self.a = {}
        self.b = {}
        
        self.a['normal'] = X[:, 0]
        self.a['tumour'] = X[:, 2]
        
        self.b['normal'] = X[:, 1]
        self.b['tumour'] = X[:, 3]
        
        self.nrows = X.shape[0]
        
    def get_independent_data(self, genome):
        return IndependentData(self.a[genome], self.b[genome])
        

class IndependentData(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b
        
        self.nrows = a.shape[0]
    
class SampleTypeException(Exception):
    pass
