'''
Define interfaces for Classifier, ClassifierRefIterator and ClassifierRow objects.

Work with convention that all C level access of position is 0-based and Python level is 1-based.

Created on 2011-07-27

@author: Andrew Roth
'''
cdef class Classifier(object):
    '''
    Base clase for all objects for classifying somatics from Bam files.
    
    Serves two functions. 
    1) Keeps track of available references accessible by refs property. 
    2) Provides iterators over references via iter_ref() method.
    '''
    def __init__(self, Counter counter):        
        self._counter = counter            
        self._refs = counter.refs
            
    def iter_ref(self, ref, **kwargs):
        '''
        Returns a ClassifierIter subclass over ref.
        '''
        raise NotImplemented
    
    property refs:
        '''
        Read only access to list of available references.
        '''
        def __get__(self):
            return self._refs

cdef class ClassifierRefIterator(object):
    '''
    Base class for all iterator objects over a reference. Should return a ClassifierRow subclass object on each
    iteration.
    '''
    def __init__(self, char * ref, JointRefIterator iter, **kwargs):
        self._ref = ref        
        self._iter = iter
        
        self._min_normal_depth = kwargs.get('min_normal_depth', 10)
        self._min_tumour_depth = kwargs.get('min_tumour_depth', 10)
    
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
        C level next method.
        
        All sub-classes need to re-implement the _get_labels() method for this to work.
        '''
        cdef int normal_depth, tumour_depth
        cdef PairedSampleBinomialCounterRow counter_row
        cdef tuple labels     
        
        while True:
            self._iter.cnext()
            
            counter_row = self._iter._current_row
            
            normal_depth = counter_row._normal_depth
            tumour_depth = counter_row._tumour_depth
            
            if normal_depth >= self._min_normal_depth and tumour_depth >= self._min_tumour_depth:
                break
        
        labels = self._get_labels()
        
        self._current_row = makeClassifierRow(counter_row, labels)
    
    cdef tuple _get_labels(self):
        pass
    
    property ref:
        '''
        Read only access to reference which the iterator runs over.
        '''
        def __get__(self):
            return self._ref

cdef class ClassifierRow(object):
    '''
    Base class for all counts row objects.
    '''
    def __str__(self):
        '''
        Method to display row object. Outputs in format tab delmited format
        
        ref position counts
        '''
        out_row = [self.ref, str(self.position), self._ref_base, self._non_ref_base]
        out_row.extend([str(x) for x in self.counts])
        out_row.extend([str(x) for x in self.labels])
        
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
            return self._counts
    
    property labels:
        '''
        Labels soft or hard over the 9 joint genotype states.
        '''
        def __get__(self):
            return self._labels
        
cdef inline ClassifierRow makeClassifierRow(PairedSampleBinomialCounterRow counter_row, tuple labels):
    '''
    Constructor method for creating a ClassifierRow from C.
    '''
    cdef ClassifierRow row = ClassifierRow.__new__(ClassifierRow)
    
    row._ref = counter_row._ref
    
    row._position = counter_row._position
    
    row._ref_base = counter_row._ref_base
    row._non_ref_base = counter_row._non_ref_base
        
    row._counts = counter_row.counts
    
    row._labels = labels
    
    return row   
