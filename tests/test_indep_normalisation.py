'''
Created on 2011-04-18

@author: andrew
'''
import unittest

import numpy as np

from joint_snv_mix.classification.indep import PairedIndependentModel
from joint_snv_mix.classification.indep_binomial import IndependentBinomialModel, IndependentBinomialParameterParser
from joint_snv_mix.classification.utils.data import JointData


class Test(unittest.TestCase):
    def setUp(self):
        param_parser = IndependentBinomialParameterParser()
        param_parser.load_from_file('../config/indep_bin.params.cfg')
        self.params = param_parser.to_dict()
    
    # Test issue 5
    def test_independent_responsibilities_normalised(self):       
        counts = np.random.randint(0, 100, size=(100, 4))
        
        data = JointData(counts)
        
        model = PairedIndependentModel(IndependentBinomialModel())
        
        resp = model.classify(data, self.params)
        
        sum_across_columns = resp.sum(axis=1)

        is_close_to_one = np.allclose(sum_across_columns, 1, 0, 1e-12)

        self.assertTrue(is_close_to_one)

        
        
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
