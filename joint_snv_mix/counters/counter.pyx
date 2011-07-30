'''
Define interfaces for Counter, CounterIter and CounterRow objects.

Work with convention that all C level access of position is 0-based and Python level is 1-based.

Created on 2011-07-06

@author: Andrew Roth
'''
cdef class Counter(object):
    '''
    Base clase for all objects for counting Bam files.
    
    Serves two functions. 
    1) Keeps track of available references accessible by refs property. 
    2) Provides iterators over references via iter_ref() method.
    '''
    def iter_ref(self, ref):
        '''
        Returns a CounterIter subclass over ref.
        '''
        raise NotImplemented
    
    property refs:
        '''
        Read only access to list of available references.
        '''
        def __get__(self):
            return self._refs
            