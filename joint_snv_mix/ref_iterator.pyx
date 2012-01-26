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
