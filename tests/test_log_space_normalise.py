'''
Created on 2010-08-31

@author: bcca
'''
import unittest

import numpy as np

from joint_snv_mix.classification.utils.normalise import log_space_normalise_rows

class Test(unittest.TestCase):
    def test_normalises_across_columns(self):
        nrows = 100
        ncols = 9
        shape = (nrows, ncols)
        log_likelihoods = np.random.random(shape)
        resp = log_space_normalise_rows(log_likelihoods)

        sum_across_columns = resp.sum(axis=1)

        is_close_to_one = np.allclose(sum_across_columns, 1, 0, 1e-12)

        self.assertTrue(is_close_to_one)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
