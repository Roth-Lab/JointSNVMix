#=======================================================================================================================
# SNVMix1 Code
#=======================================================================================================================
cdef double * multinomial_log_likelihood(int * counts,
                                         double ** log_mu,
                                         int num_genotypes,
                                         int num_bases):
    '''
    Return class log_likelihoods under SNVMix1 model with parameters mu=exp(log_mu).
    
    Allocates log_likelihood which will need to be freed by caller.
    '''
    cdef int b, g
    cdef double * log_likelihood = < double *> malloc(num_genotypes * sizeof(double))
    
    for g in range(num_genotypes):
        log_likelihood[g] = 0
        
        for b in range(num_bases):
            log_likelihood[g] += counts[b] * log_mu[g][b] 
    
    return log_likelihood

#=======================================================================================================================
# SNVMix2 Code
#=======================================================================================================================
cdef double * get_phred_to_prob_qual_map(int num_quals):
        cdef int qual
        cdef double base, exp
        cdef double * qual_map = < double *> malloc(num_quals * sizeof(double))
        
        for qual in range(num_quals):
            exp = -1 * (< double > qual) / 10
            base = 10
            
            qual_map[qual] = 1 - pow(base, exp)
        
        return qual_map

cdef double * snv_mix_2_log_likelihood(base_map_qualities_struct data,
                                       double * mu,
                                       double * qual_map,
                                       int num_genotypes):
    '''
    Return class log_likelihoods under SNVMix2 model with parameters mu.
    
    Allocates log_likelihood which will need to be freed by caller.
    '''
    cdef int genotype, read_index
    cdef double ll, m
    cdef double * log_likelihood = < double *> malloc(num_genotypes * sizeof(double)) 
    
    for genotype in range(num_genotypes):
        ll = 0            
        m = mu[genotype]
        
        for read_index in range(data.depth.A):
             ll += compute_single_match_base_log_prob(data.base_quals.A[read_index],
                                                      data.map_quals.A[read_index],
                                                      m,
                                                      qual_map)
        
        for read_index in range(data.depth.B):
             ll += compute_single_mismatch_base_log_prob(data.base_quals.B[read_index],
                                                         data.map_quals.B[read_index],
                                                         m,
                                                         qual_map)
        
        log_likelihood[genotype] = ll
    
    return log_likelihood
        

cdef double compute_single_match_base_log_prob(int base_qual, int map_qual, double mu, double * qual_map):
    cdef double q, r   

    q = qual_map[base_qual]
    r = qual_map[map_qual]
    
    return compute_single_base_log_prob(q, r, mu)         

cdef double compute_single_mismatch_base_log_prob(int base_qual, int map_qual, double mu, double * qual_map):
    cdef double q, r    

    q = (1 - qual_map[base_qual]) / 3
    r = qual_map[map_qual]
    
    return compute_single_base_log_prob(q, r, mu)    
        
cdef double compute_single_base_log_prob(double q, double r, double m):
    aligned_wrong = 0.5 * (1 - r)
    aligned_right = r * ((1 - q) * (1 - m) + q * m)
    
    return log(aligned_wrong + aligned_right)

#=======================================================================================================================
# Code for joining independent probabilities to joint probabilities.
#=======================================================================================================================
cdef double * combine_independent_probs(double * normal_probs,
                                        double * tumour_probs,
                                        int num_normal_genotypes,
                                        int num_tumour_genotypes):
    '''
    Combine independent probabilities from a normal and tumour sample to get joint probabilities.
    
    Allocates joint_probabilities which needs to be freed by caller.
    '''    
    cdef int normal_index, tumour_index, joint_index, num_joint_genotypes
    cdef double total
    cdef double * joint_probs
    
    num_joint_genotypes = num_normal_genotypes * num_tumour_genotypes
    joint_probs = < double *> malloc(num_joint_genotypes * sizeof(double))
    
    total = 0
    
    for normal_index in range(num_normal_genotypes):
        for tumour_index in range(num_tumour_genotypes):
            joint_index = num_normal_genotypes * normal_index + tumour_index
            
            joint_probs[joint_index] = normal_probs[normal_index] * tumour_probs[tumour_index]

            total += joint_probs[joint_index]
            
    for joint_index in range(num_joint_genotypes):
        joint_probs[joint_index] = joint_probs[joint_index] / total
    
    return joint_probs

#=======================================================================================================================
# General mixture model code.
#=======================================================================================================================
cdef double * get_mixture_posterior(double * log_likelihood,
                                    double * log_mix_weight,
                                    int num_classes):
    '''
    Compute normalised posterior probabilities from mixture distribution by adding log_likelihood and 
    then normalising.
    
    Allocates posterior which will need to be freed by caller.
    ''' 
    cdef int i
    cdef double * posterior = < double *> malloc(num_classes * sizeof(double))

    for i in range(num_classes):
        log_likelihood[i] = log_mix_weight[i] + log_likelihood[i]
    
    log_space_normalise_row(log_likelihood, num_classes)
    
    for i in range(num_classes):
        posterior[i] = exp(log_likelihood[i])
    
    return posterior

cdef double * get_joint_posterior(double * normal_log_likelihood,
                                  double * tumour_log_likelihood,
                                  double * log_mix_weight,
                                  int num_normal_genotypes,
                                  int num_tumour_genotypes):
    
    cdef int normal_index, tumour_index, joint_index, num_joint_genotypes
    cdef double total
    cdef double * joint_log_likelihood, * joint_posterior
    
    num_joint_genotypes = num_normal_genotypes * num_tumour_genotypes
    joint_log_likelihood = < double *> malloc(num_joint_genotypes * sizeof(double)) 
    
    for normal_index in range(num_normal_genotypes):
        for tumour_index in range(num_tumour_genotypes):
            joint_index = num_normal_genotypes * normal_index + tumour_index
            
            joint_log_likelihood[joint_index] = normal_log_likelihood[normal_index] + tumour_log_likelihood[tumour_index]
            
    joint_posterior = get_mixture_posterior(joint_log_likelihood,
                                            log_mix_weight,
                                            num_joint_genotypes)
    
    free(joint_log_likelihood)
    
    return joint_posterior

#=======================================================================================================================
# Code for doing log space normalisation
#=======================================================================================================================
cdef void log_space_normalise_row(double * log_X, int size):
    '''
    Normalise log_X so that 
    
    exp(log_X[0]) + ... + exp(log_X[1]) == 1
    
    Done in place so log_X is modified.
    '''
    cdef int i
    cdef double norm_const
    
    norm_const = log_sum_exp(log_X, size)
    
    for i in range(size):
        log_X[i] = log_X[i] - norm_const    

cdef double log_sum_exp(double * log_X, int size):
    '''
    Given a c-array log_X of values compute log( exp(log_X[0]) + ... + exp(log_X[size]) ).
    
    Numerically safer than naive method.
    '''
    cdef int i
    cdef double max_exp, total
 
    max_exp = log_X[0]
 
    for i in range(size):
        if max_exp < log_X[i]:
            max_exp = log_X[i]

    total = 0
    for i in range(size):
        total += exp(log_X[i] - max_exp)
    
    return log(total) + max_exp
