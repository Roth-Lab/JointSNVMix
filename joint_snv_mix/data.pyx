'''
Created on 2012-01-16

@author: Andrew Roth
'''
#=======================================================================================================================
# Data-sets
#=======================================================================================================================
class DataSet(object):
    def __init__(self):
        self._data_points = []
    
    def __iter__(self):
        return iter(self._data_points)
    
    def add_point(self, data_point):
        self._data_points.append(data_point)    

#=======================================================================================================================
# Data-points
#=======================================================================================================================
cdef class DataPoint(object):
    '''
    Base class for all data point objects.
    '''
    pass

cdef class BiAllelicCountDataPoint(object):
    '''
    Class for representing A/B ref/non-ref allelic count data.
    '''
    pass
        
cdef class BiAllelicQualityDataPoint(object):
    '''    
    Class for representing A/B ref/non-ref allelic alignment and mapping quality data.
    '''
    def __dealloc__(self):
        free(self._q_N)
        free(self._q_T)
        free(self._r_N)
        free(self._r_T)

#=======================================================================================================================
# Data-point factory functions
#=======================================================================================================================
cdef BiAllelicCountDataPoint makeBiAllelicCountDataPoint(int a_N, int b_N, int a_T, int b_T):
    cdef BiAllelicCountDataPoint data_point = BiAllelicCountDataPoint.__new__(BiAllelicCountDataPoint)
    
    data_point._a_N = a_N
    data_point._b_N = b_N
    
    data_point._a_T = a_T
    data_point._b_T = b_T
    
    return data_point

cdef BiAllelicQualityDataPoint makeBiAllelicQualityDataPoint(int num_normal_reads, int num_tumour_reads,
                                                             double * q_N, double * r_N,
                                                             double * q_T, double * r_T):
    cdef BiAllelicQualityDataPoint data_point = BiAllelicQualityDataPoint.__new__(BiAllelicQualityDataPoint)
    
    data_point._num_normal_reads = num_normal_reads
    data_point._num_tumour_reads = num_tumour_reads
    
    data_point._q_N = q_N
    data_point._r_N = r_N
    
    data_point._q_T = q_T
    data_point._r_T = b_T
    
    return data_point
