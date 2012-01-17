'''
Created on 2012-01-16

@author: Andrew Roth
'''
cdef class Model(object):
    '''
    Base class which all models inherit from.
    '''
    def __init__(self, *args, **kwargs):
        pass
    
    cdef fit(self, Dataset data, **kwargs):
        '''
        Function to fit the model to the dataset. Should return a valid Parameter object.
        '''
        pass
    
    cdef predict(self, DataPoint data, **kwargs):
        '''
        Function to predict the class which the data point belongs to. Should return a Prediction object.
        '''
        pass
    
