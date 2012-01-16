#cython: cdivision=True

DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class SnvMixTwoClassifier(Classifier):
    def __init__(self, **kwargs):       
        self._params = PairedSnvMixParameters(**kwargs)

    cdef double * _get_labels(self, PairedSampleBinomialCounterRow row):
        cdef SnvMixTwoData normal_data, tumour_data
        cdef SnvMixTwoCpt normal_cpt, tumour_cpt
        cdef double * normal_resp, * tumour_resp, * joint_resp
        
        normal_data = makeSnvMixTwoData((< JointBinaryQualityCounterRow > row)._normal_data)
        tumour_data = makeSnvMixTwoData((< JointBinaryQualityCounterRow > row)._tumour_data)
        
        normal_cpt = SnvMixTwoCpt(normal_data, self._params._normal_params)
        tumour_cpt = SnvMixTwoCpt(tumour_data, self._params._tumour_params)
        
        normal_resp = normal_cpt.get_resp()
        tumour_resp = tumour_cpt.get_resp()
                        
        joint_resp = combine_independent_probs(normal_resp,
                                               tumour_resp,
                                               NUM_GENOTYPES,
                                               NUM_GENOTYPES)        
        free(normal_resp)
        free(tumour_resp)
               
        return joint_resp
        
