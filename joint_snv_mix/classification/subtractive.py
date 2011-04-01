'''
Created on 2011-03-31

@author: andrew
'''
def get_joint_log_likelihoods(log_likelihoods, pi):
    normal_log_likelihoods = log_likelihoods['normal']
    tumour_log_likelihoods = log_likelihoods['tumour']
    
    normal_nclass = normal_log_likelihoods.shape[1]
    column_shape = (normal_log_likelihoods[:, 0].size, 1)

    log_likelihoods = np.hstack([normal_log_likelihoods[:, i].reshape(column_shape) + tumour_log_likelihoods
                                  for i in range(normal_nclass)])
    
    log_likelihoods = log_likelihoods + np.log(pi)

    return log_likelihoods