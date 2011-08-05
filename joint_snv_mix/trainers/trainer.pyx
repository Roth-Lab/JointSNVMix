'''
Created on 2011-08-04

@author: Andrew Roth
'''
cdef class Trainer(object):   
    cdef list _load_data_set(self, Counter counter, double sub_sample_fraction):
        cdef int ref_len
        cdef RefIterator ref_iter
        cdef gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937)
        cdef list sample_data = []
        cdef PairedSampleBinomialCounterRow row
        
        print "Randomly sub-sampling {0}% of the data set.".format(100 * sub_sample_fraction)
        
        for ref in ['1', ]:
            ref_iter = counter.iter_ref(ref)
            
            ref_len = 0
            
            try:
                while True:
                    ref_iter.cnext()
                    
                    row = ref_iter._current_row
                    
                    if row._normal_depth < 10 or row._tumour_depth < 10:
                        continue
                    
                    if gsl_rng_uniform(r) <= sub_sample_fraction:                        
                        sample_data.append(row)
                        
                        ref_len += 1
                        
            except StopIteration:
                pass
            
            print "Sub-sampled {0} positions from ref {1}".format(ref_len, ref)
        
        print "Total sub-sample size is {0}".format(len(sample_data))
        
        return sample_data
