'''
Created on 2012-01-25

@author: Andrew Roth
'''
cdef class RefIterator(object):
    def __iter__(self):
        return self
    
    def __next__(self):
        self.cnext()
        
        return self._current_row

    cdef cnext(self):
        self.advance_position()
        self.parse_current_position()
        
    cdef advance_position(self):
        pass
    
    cdef parse_current_position(self):
        pass
    
    property ref:
        '''
        Read only access to reference which the iterator runs over.
        '''
        def __get__(self):
            return self._ref
    
    property position:
        '''
        Read only access to 1-based current position of iterator.
        '''
        def __get__(self):
            return self._pos + 1     
