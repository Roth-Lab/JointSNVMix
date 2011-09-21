class SubSampler(object):    
    def __init__(self, skip_size, min_normal_depth, min_tumour_depth):
        self._sampler = PairedDataSubSampler(skip_size, min_normal_depth, min_tumour_depth)

class SnvMixOneSubSampler(SubSampler):        
    def subsample(self, counter, refs=None):
        sample = SnvMixOneSubsample()
        
        self._sampler.subsample(counter, sample, refs)
        
class SnvMixTwoSubSampler(SubSampler):        
    def subsample(self, counter, refs=None):
        sample = SnvMixTwoSubsample()
        
        self._sampler.subsample(counter, sample, refs)        
                
#---------------------------------------------------------------------------------------------------------------------- 
cdef class PairedDataSubSampler(object): 
    def __init__(self, int skip_size, int min_normal_depth, int min_tumour_depth):
        self._skip_size = skip_size
        self._min_normal_depth = min_normal_depth
        self._min_tumour_depth = min_tumour_depth

    def subsample(self, Counter counter, Subsample sample, refs=None):
        cdef int i, ref_sample_size
        cdef RefIterator ref_iter       
        cdef PairedSampleBinomialCounterRow row
        
        print '''Randomly sub-sampling every {0}th position with 
                 normal depth {1} and tumour depth {2} the data set.'''.format(self._skip_size,
                                                                               self._min_normal_depth,
                                                                               self._min_tumour_depth)
        
        if refs == None:
            refs = counter.refs
        
        for ref in refs:
            print "Subsampling ref {0}.".format(ref)

            ref_iter = counter.iter_ref(ref)
            
            i = 0
            ref_sample_size = 0
            
            try:
                while True:
                    ref_iter.cnext()
                    
                    row = ref_iter._current_row
                    
                    if row._normal_depth < self._min_normal_depth or row._tumour_depth < self._min_tumour_depth:
                        continue
                    
                    if i % self._skip_size == 0:               
                        sample.add_row(row)
                        ref_sample_size += 1
                        
                    i += 1
                        
            except StopIteration:
                pass
            
            print "Sub-sampled {0} positions from ref {1}".format(ref_sample_size, ref)
        
        print "Total sub-sample size is {0}".format(i)
        
        return sample

#---------------------------------------------------------------------------------------------------------------------- 
cdef class Subsample(object):
    def __init__(self):
        self._normal_rows = []
        self._tumour_rows = []
    
    cdef add_row(self, PairedSampleBinomialCounterRow row):
        pass
#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixOneSubsample(Subsample):
    cdef add_row(self, PairedSampleBinomialCounterRow row):
        cdef SnvMixOneData normal_data, tumour_data
    
        normal_data = makeSnvMixOneData((< JointBinaryCounterRow > row)._normal_counts)
        tumour_data = makeSnvMixOneData((< JointBinaryCounterRow > row)._tumour_counts)
        
        self._normal_rows.append(normal_data)
        self._tumour_rows.append(tumour_data)

#---------------------------------------------------------------------------------------------------------------------- 
cdef class SnvMixTwoSubsample(Subsample):
    cdef add_row(self, PairedSampleBinomialCounterRow row):
        cdef SnvMixTwoData normal_data, tumour_data
    
        normal_data = makeSnvMixTwoData((< JointBinaryQualityCounterRow > row)._normal_data)
        tumour_data = makeSnvMixTwoData((< JointBinaryQualityCounterRow > row)._tumour_data)
        
        self._normal_rows.append(normal_data)
        self._tumour_rows.append(tumour_data)
