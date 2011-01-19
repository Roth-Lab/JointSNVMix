'''
Created on 2011-01-18

@author: Andrew Roth
'''
#=======================================================================================================================
# Independent Models
#=======================================================================================================================
def independent_binomial_log_likelihood( data, parameters ):
    a = data.a
    b = data.b
    
    d = a + b
    
    mu = parameters['mu']

    log_likelihoods = log_binomial_likelihood( a, d, mu )

    pi = parameters['pi']
    log_pi = np.log( pi )

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods

def independent_binomial_log_likelihood( data, parameters ):
    a = data.a
    b = data.b
    
    d = a + b
    
    mu = parameters['mu']

    log_likelihoods = log_binomial_likelihood( a, d, mu )

    pi = parameters['pi']
    log_pi = np.log( pi )

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods

#=======================================================================================================================
# Joint Models
#=======================================================================================================================
def joint_beta_binomial_log_likelihood( data, parameters ):
    a_1 = data.a[0]
    a_2 = data.a[1]
    
    d_1 = data.a[0] + data.b[0]
    d_2 = data.a[1] + data.b[1]
    
    alpha_1 = parameters['alpha'][0]
    alpha_2 = parameters['alpha'][1]
    
    beta_1 = parameters['beta'][0]
    beta_2 = parameters['beta'][1]
    
    normal_log_likelihoods = log_beta_binomial_likelihood( a_1, d_1, alpha_1, beta_1 )
    tumour_log_likelihoods = log_beta_binomial_likelihood( a_2, d_2, alpha_2, beta_2 )

    column_shape = ( normal_log_likelihoods[:, 0].size, 1 )

    log_likelihoods = np.hstack( ( 
                                 normal_log_likelihoods[:, 0].reshape( column_shape ) + tumour_log_likelihoods ,
                                 normal_log_likelihoods[:, 1].reshape( column_shape ) + tumour_log_likelihoods ,
                                 normal_log_likelihoods[:, 2].reshape( column_shape ) + tumour_log_likelihoods
                                 ) )

    pi = parameters['pi']
    log_likelihoods = log_likelihoods + np.log( pi )
    
    return log_likelihoods

def joint_binomial_log_likelihood( data, parameters ):
    a_1 = data.a[0]
    a_2 = data.a[1]
    
    d_1 = data.a[0] + data.b[0]
    d_2 = data.a[1] + data.b[1]
    
    mu_1 = parameters['mu'][0]
    mu_2 = parameters['mu'][1]

    
    normal_log_likelihoods = log_binomial_likelihood( a_1, d_1, mu_1 )
    tumour_log_likelihoods = log_binomial_likelihood( a_2, d_2, mu_2 )

    column_shape = ( normal_log_likelihoods[:, 0].size, 1 )

    log_likelihoods = np.hstack( ( 
                                 normal_log_likelihoods[:, 0].reshape( column_shape ) + tumour_log_likelihoods ,
                                 normal_log_likelihoods[:, 1].reshape( column_shape ) + tumour_log_likelihoods ,
                                 normal_log_likelihoods[:, 2].reshape( column_shape ) + tumour_log_likelihoods
                                 ) )

    pi = parameters['pi']
    log_pi = np.log( pi )

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods
