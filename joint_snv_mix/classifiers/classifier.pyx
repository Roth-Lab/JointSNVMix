'''
Define interfaces for Classifier and ClassifierRow objects.

Created on 2011-07-27

@author: Andrew Roth
'''
cdef class Classifier(object):
    def classify(self, row):
        return self._classify(row)
        
    cdef ClassifierRow _classify(self, PairedSampleBinomialCounterRow row):
        cdef tuple labels
        
        labels = self._get_labels(row)
        
        return makeClassifierRow(row, labels)
    
    cdef tuple _get_labels(self, PairedSampleBinomialCounterRow row):
        pass

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
        out_row.extend(["{0:.4f}".format(x) for x in self.labels])
        
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
