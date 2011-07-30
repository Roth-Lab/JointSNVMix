'''
Created on 2011-07-29

@author: Andrew Roth
'''
cdef class CounterRow(object):
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")
    
    property ref:
        def __get__(self):
            return self._ref
    
    property position:
        '''
        1-based position
        '''
        def __get__(self):
            return self._position + 1

cdef class SingleSampleCounterRow(CounterRow):    
    property depth:
        def __get__(self):    
            return self._depth

cdef class PairedSampleCounterRow(object):
    property normal_depth:
        def __get__(self):
            return self._normal_depth
    
    property tumour_depth:
        def __get__(self):
            return self._tumour_depth