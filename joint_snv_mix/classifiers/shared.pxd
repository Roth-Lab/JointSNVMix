from libc.math cimport exp, log
from libc.stdlib cimport malloc, free

from joint_snv_mix.counters.joint_binary_quality_counter cimport base_map_qualities_struct

cdef double multinomial_log_likelihood(int * counts, double * log_mu, int num_bases)

cdef double * get_phred_to_prob_qual_map(int num_quals)

cdef double * snv_mix_2_log_likelihood(base_map_qualities_struct data,
                                       double * mu,
                                       double * qual_map,
                                       int num_genotypes)

cdef double * combine_independent_probs(double * normal_probs,
                                        double * tumour_probs,
                                        int num_normal_genotypes,
                                        int num_tumour_genotypes)

cdef double * get_mixture_posterior(double * log_likelihood, double * log_mix_weight, int num_classes)

cdef double * get_joint_posterior(double * normal_log_likelihood,
                                  double * tumour_log_likelihood,
                                  double * log_mix_weight,
                                  int num_normal_genotypes,
                                  int num_tumour_genotypes)
