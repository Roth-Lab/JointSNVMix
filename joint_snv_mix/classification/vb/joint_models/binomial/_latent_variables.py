'''
Created on 2010-11-17

@author: Andrew Roth
'''
import numpy as np

class _LatentVariables(object):
    def __init__(self, data,):
        nrows = data.shape[0]
        shape = (nrows, 1)
        
        self.a_1 = data[:, 0].reshape(shape)
        self.a_2 = data[:, 2].reshape(shape)
        
        self.b_1 = data[:, 1] - data[:, 0]
        self.b_1 = self.b_1.reshape(shape)
        
        self.b_2 = data[:, 3] - data[:, 2]
        self.b_2 = self.b_2.reshape(shape)

    def update(self, Posterior):
        term_1 = self.a_1 * Posterior.expected_log_mu_1 + \
                 self.b_1 * Posterior.expected_log_one_minus_mu_1
        
        term_2 = self.a_2 * Posterior.expected_log_mu_2 + \
                 self.b_2 * Posterior.expected_log_one_minus_mu_2

        self.responsibilities = np.column_stack((
                                                 Posterior.expected_log_pi[0] + term_1[:, 0] + term_2[:, 0],
                                                 Posterior.expected_log_pi[1] + term_1[:, 0] + term_2[:, 1],
                                                 Posterior.expected_log_pi[2] + term_1[:, 0] + term_2[:, 2],
                                                 Posterior.expected_log_pi[3] + term_1[:, 1] + term_2[:, 0],
                                                 Posterior.expected_log_pi[4] + term_1[:, 1] + term_2[:, 1],
                                                 Posterior.expected_log_pi[5] + term_1[:, 1] + term_2[:, 2],
                                                 Posterior.expected_log_pi[6] + term_1[:, 2] + term_2[:, 0],
                                                 Posterior.expected_log_pi[7] + term_1[:, 2] + term_2[:, 1],
                                                 Posterior.expected_log_pi[8] + term_1[:, 2] + term_2[:, 2]
                                                 ))

        self._normalise()
        
        self.marginal_responsibilities_1 = np.column_stack((
                                                          self.responsibilities[:, 0] + self.responsibilities[:, 1] + self.responsibilities[:, 2],
                                                          self.responsibilities[:, 3] + self.responsibilities[:, 4] + self.responsibilities[:, 5],
                                                          self.responsibilities[:, 6] + self.responsibilities[:, 7] + self.responsibilities[:, 8]
                                                          ))

        self.marginal_responsibilities_2 = np.column_stack((
                                                          self.responsibilities[:, 0] + self.responsibilities[:, 3] + self.responsibilities[:, 6],
                                                          self.responsibilities[:, 1] + self.responsibilities[:, 4] + self.responsibilities[:, 7],
                                                          self.responsibilities[:, 2] + self.responsibilities[:, 5] + self.responsibilities[:, 8]
                                                          ))
                                

    def _normalise(self):
        nrows = self.responsibilities.shape[0]
        
        log_norm_const = np.logaddexp.reduce(self.responsibilities, axis=1)
        log_norm_const = log_norm_const.reshape(nrows, 1)

        self.responsibilities = self.responsibilities - log_norm_const
        self.responsibilities = np.exp(self.responsibilities)
        
        eps = np.finfo(np.float64).eps
        
        self.responsibilities[self.responsibilities <= eps] = 0



