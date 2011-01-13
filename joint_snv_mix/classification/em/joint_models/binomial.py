'''
Created on 2010-12-09

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix.classification.em.em_lower_bound import EMLowerBound
from joint_snv_mix.classification.em.em_model import EMModel, EMModelTrainer
from joint_snv_mix.classification.em.em_posterior import EMPosterior
from joint_snv_mix.classification.em.joint_models.joint_latent_variables import JointLatentVariables
from joint_snv_mix.classification.utils.log_pdf import log_beta_pdf, log_binomial_likelihood
from joint_snv_mix.classification.utils.marginal_responsibilities import get_marginals

class JointBinomialModel( EMModel ):
    def __init__( self ):
        self.trainer_class = JointBinomialModelTrainer
        
        self.log_likelihood_func = joint_binomial_log_likelihood

class JointBinomialModelTrainer( EMModelTrainer ):
    def _init_components( self ):
        self.latent_variables = JointBinomialLatentVariables( self.data )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = JointBinomialPosterior( self.data, self.priors, self.responsibilities )
        
        self.lower_bound = JointBinomialLowerBound( self.data, self.priors )

class JointBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_binomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        log_prior = 0.
        
        for genome in range( 2 ):
            mu = self.parameters['mu'][genome]
            alpha = self.priors['alpha'][genome]
            beta = self.priors['beta'][genome]
            
            log_prior += log_beta_pdf( mu, alpha, beta )
        
        log_prior = log_prior.sum()

        return log_prior
    
class JointBinomialLatentVariables( JointLatentVariables ):
    def __init__( self, data ):
        JointLatentVariables.__init__( self, data )
        
        self.likelihood_func = joint_binomial_log_likelihood
        
class JointBinomialPosterior( EMPosterior ):
    def _init_parameters( self ):
        self.parameters = {}
        
        self._update_mix_weights()
        
        self.parameters['mu'] = []
        
        for i in range( 2 ):
            alpha = self.priors['alpha'][i, :]
            beta = self.priors['beta'][i, :]
            
            mu = alpha / ( alpha + beta )
            
            self.parameters['mu'].append( mu )
           
    def _update_density_parameters( self ):
        marginals = get_marginals( self.responsibilities )
        
        for genome in range( 2 ):
            a = self.data.a[genome]
            b = self.data.b[genome]
            
            alpha = self.priors['alpha'][genome]
            beta = self.priors['beta'][genome]
            
            tau = marginals[genome]
            
            self.parameters['mu'][genome] = self._update_mu( a, b, alpha, beta, tau )
    
    def _update_mu( self, a, b, alpha, beta, tau ):       
        d = a + b
        
        n = self.data.nrows
        shape = ( n, 1 )
        
        a = a.reshape( shape )
        d = d.reshape( shape )
        
        ref_sum = np.sum( tau * a, axis=0 )
        
        depth_sum = np.sum( tau * d, axis=0 )
        
        numerator = ref_sum + alpha - 1.
        
        denominator = depth_sum + alpha + beta - 2.
        
        return np.exp( np.log( numerator ) - np.log( denominator ) )
        
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
