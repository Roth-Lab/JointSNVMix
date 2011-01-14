'''
Created on 2010-11-17

@author: Andrew Roth
'''
import numpy as np

class _ESS(object):
    def __init__(self, data):
        nrows = data.nrows
        shape = ( nrows, 1 )
        
        self.a_1 = data.a[0].reshape( shape )
        self.a_2 = data.a[1].reshape( shape )
        
        self.b_1 = data.b[0].reshape( shape )
        self.b_2 = data.b[1].reshape( shape )
        
    def update(self, LatentVariables):
        self.N_g = LatentVariables.responsibilities.sum(axis=0)
        
        self.a_bar_1 = LatentVariables.marginal_responsibilities_1 * self.a_1
        self.a_bar_1 = self.a_bar_1.sum(axis=0)
        
        self.a_bar_2 = LatentVariables.marginal_responsibilities_2 * self.a_2
        self.a_bar_2 = self.a_bar_2.sum(axis=0)
        
        self.b_bar_1 = LatentVariables.marginal_responsibilities_1 * self.b_1
        self.b_bar_1 = self.b_bar_1.sum(axis=0)
        
        self.b_bar_2 = LatentVariables.marginal_responsibilities_2 * self.b_2
        self.b_bar_2 = self.b_bar_2.sum(axis=0)
