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

cdef class CounterRefIterator(object):
    '''
    Base class for all iterator objects over a reference. Should return a CounterRow subclass object on each iteration.
    '''
    def __iter__(self):
        return self
    
    def __next__(self):
        '''
        Python level next() method.
        '''
        self.cnext()
        
        return self._current_row
    
    cdef cnext(self):
        '''
        C level next method. All subclasses need to re-implement this method.
        
        Should move the self._current_row to next value.
        
        No return value should be specified to allow Cython to pass exceptions along chain.
        '''
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
            return self._position + 1

cdef class CounterRow(object):
    '''
    Base class for all counts row objects.
    '''
    def __str__(self):
        '''
        Method to display row object. Outputs in format tab delmited format
        
        ref position counts
        '''
        out_row = [self.ref, str(self.position)]
        out_row.extend([str(x) for x in self.counts])
        
        return "\t".join(out_row)
    
    property ref:
        '''
        Read only access to reference for row.
        '''
        def __get__(self):
            return self._ref
    
    property position:
        '''
        Read only access to 1-based position of row.
        '''
        def __get__(self):
            return self._position + 1
    
    property counts:
        '''
        Read only access to count data. Should return the counts as a list of integers.
        '''
        def __get__(self):
            pass
    
    property depth:
        '''
        Read only access to depth. Depth will depend on the nature of CounterRow. For rows with counts from multiple
        references minimum depth over samples should be returned.
        '''
        def __get__(self):
            pass
            