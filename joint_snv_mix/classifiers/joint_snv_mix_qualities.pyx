DEF NUM_GENOTYPES = 3
DEF NUM_JOINT_GENOTYPES = 9

cdef class JointSnvMixTwoClassifier(Classifier):
    def __init__(self, **kwargs):
        self._params = JointSnvMixParameters(**kwargs)

    cdef double * _get_labels(self, PairedSampleBinomialCounterRow row):
        cdef JointSnvMixTwoCpt cpt
        cdef JointSnvMixTwoData data
        
        data = makeJointSnvMixTwoData(row) 
        
        cpt = JointSnvMixTwoCpt(data, self._params)
        
        return cpt.get_resp()
