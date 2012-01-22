'''
Created on 2012-01-22

@author: Andrew Roth
'''
        
results_header = ['chrom',
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
                  'p_BB_BB']

cdef class CResultsWriter(object):
    def __init__(self, file_name=None):
        if file_name is None:
            self._file_ptr = stdout
        else:
            self._file_ptr = fopen(< char *> file_name, "w")
        
            if self._file_ptr == NULL:
                raise Exception("Couldn't open results file at {0}".format(file_name))
        
        self._init_format_string()
        
        header = "\t".join(results_header)
        header += "\n"
        
        fputs(< char *> header, self._file_ptr)
    
    def __dealloc__(self):
        if self._file_ptr != NULL:
            self.close()
        
        if self._format_string != NULL:
            free(self._format_string)
            self._format_string = NULL
    
    cdef _init_format_string(self):
        format_string = ['%s', # chrom
                         '%d', # position
                         '%s', # ref_base
                         '%s', # var_base
                         '%d', # normal_ref_counts
                         '%d', # normal_var_counts
                         '%d', # tumour_ref_counts
                         '%d', # tumour_var_counts
                         '%.4f', # p_AA_AA
                         '%.4f', # p_AA_AB
                         '%.4f', # p_AA_BB
                         '%.4f', # p_AB_AA
                         '%.4f', # p_AB_AB
                         '%.4f', # p_AB_BB
                         '%.4f', # p_BB_AA
                         '%.4f', # p_BB_AB
                         '%.4f', # p_BB_BB
                         ]
        
        format_string = "\t".join(format_string) + "\n"
        
        self._format_string = strdup(< char *> format_string)
        
    cdef close(self):
        fclose(self._file_ptr)
        self._file_ptr = NULL

    cdef write_position(self, JointBinaryCounterRow row, double * probs):
        cdef JointBinaryData data
        
        data = row._data
        
        fprintf(self._file_ptr, self._format_string, row._ref,
                                                     row._pos + 1,
                                                     row._ref_base,
                                                     row._var_base,
                                                     data._a_N,
                                                     data._b_N,
                                                     data._a_T,
                                                     data._b_T,
                                                     probs[0],
                                                     probs[1],
                                                     probs[2],
                                                     probs[3],
                                                     probs[4],
                                                     probs[5],
                                                     probs[6],
                                                     probs[7],
                                                     probs[8])        

