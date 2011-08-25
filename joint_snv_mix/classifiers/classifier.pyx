'''
Define interfaces for Classifier and ClassifierRow objects.

Created on 2011-07-27

@author: Andrew Roth
'''
DEF NUM_JOINT_GENOTYPES = 9

def classify_counter(Counter counter, Classifier classifier, char * file_name):
    cdef PairedSampleBinomialCounterRow row
    cdef ClassifierRow cl_row
    cdef FILE * file_p
    
    file_p = get_out_file(file_name)
    
    for ref in sorted(counter.refs):
        print ref
        
        for row in counter.iter_ref(ref):
            cl_row = classifier._classify(row)
            
            fprintf(file_p,
                    "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                    cl_row._ref,
                    cl_row._position + 1,
                    cl_row._ref_base,
                    cl_row._non_ref_base,
                    < int > cl_row._counts[0],
                    < int > cl_row._counts[1],
                    < int > cl_row._counts[2],
                    < int > cl_row._counts[3],
                    < double > cl_row._labels[0],
                    < double > cl_row._labels[1],
                    < double > cl_row._labels[2],
                    < double > cl_row._labels[3],
                    < double > cl_row._labels[4],
                    < double > cl_row._labels[5],
                    < double > cl_row._labels[6],
                    < double > cl_row._labels[7],
                    < double > cl_row._labels[8]
                    )
    
    fclose(file_p)

cdef FILE * get_out_file(char * file_name):
    cdef FILE * file_p
    
    file_p = fopen(file_name, "w")
    
    if file_p == NULL:
        raise Exception("Couldn't open out file")
    
    header = [
              'chrom',
              'position',
              'ref_base',
              'var_base',
              'normal_counts_a',
              'normal_counts_b',
              'tumour_counts_a',
              'tumour_counts_b',
              'p_AA_AA',
              'p_AA_AB',
              'p_AA_BB',
              'p_AB_AA',
              'p_AB_AB',
              'p_AB_BB',
              'p_BB_AA',
              'p_BB_AB',
              'p_BB_BB',
              '\n'
              ]
    
    header = "\t".join(header)
    fputs(< char *> header, file_p)
    
    return file_p

#=======================================================================================================================
# Objects
#=======================================================================================================================
cdef class Classifier(object):
    def classify(self, row):
        return self._classify(row)
        
    cdef ClassifierRow _classify(self, PairedSampleBinomialCounterRow row):
        cdef double * labels
        
        labels = self._get_labels(row)
        
        return makeClassifierRow(row, labels)
    
    cdef double * _get_labels(self, PairedSampleBinomialCounterRow row):
        pass

cdef class ClassifierRow(object):
    '''
    Base class for all counts row objects.
    '''
    def __dealloc__(self):
        free(self._labels)
    
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
            return [x for x in self._counts[:4]]
    
    property labels:
        '''
        Labels soft or hard over the 9 joint genotype states.
        '''
        def __get__(self):
            return [x for x in self._labels[:NUM_JOINT_GENOTYPES]]
        
cdef inline ClassifierRow makeClassifierRow(PairedSampleBinomialCounterRow counter_row, double * labels):
    '''
    Constructor method for creating a ClassifierRow from C.
    '''
    cdef int i
    cdef tuple counts
    
    cdef ClassifierRow row = ClassifierRow.__new__(ClassifierRow)
    
    row._ref = counter_row._ref
    
    row._position = counter_row._position
    
    row._ref_base = counter_row._ref_base
    row._non_ref_base = counter_row._non_ref_base
    
    for i in range(4):
        row._counts[i] = counter_row._counts[i]
    
    row._labels = labels
    
    return row   
