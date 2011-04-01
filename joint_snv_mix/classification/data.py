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

class IndependentData(object):
    def __init__(self, X, type):
        if type == 'normal':
            self.a = X[:, 0]
            self.b = X[:, 1]
        elif type == 'tumour':
            self.a = X[:, 2] 
            self.b = X[:, 3]
        else:
            raise SampleTypeException
        
        self.nrows = X.shape[0]
        
class MultinomialData(object):
    def __init__(self, X):
        self.counts = {}

        self.counts['normal'] = X[:, :4]
        self.counts['tumour'] = X[:, 4:]

        self.nrows = X.shape[0]
        
class JointQualityData(object):
    def __init__(self, X, normal_base_qualities, tumour_base_qualities):
        self.a = []
        self.b = []
        
        self.a.append(X[:, 0])
        self.a.append(X[:, 2])
        
        self.b.append(X[:, 1])
        self.b.append(X[:, 3])
        
        self.nrows = X.shape[0]
        
        self.qualities = []
        self.qualities.append(normal_base_qualities)
        self.qualities.append(tumour_base_qualities)
    
    
class SampleTypeException(Exception):
    pass
