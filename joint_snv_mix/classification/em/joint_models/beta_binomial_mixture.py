'''
Created on 2010-12-09

@author: Andrew Roth
'''
from joint_snv_mix.classification.em.em_lower_bound import EMLowerBound
from joint_snv_mix.classification.em.em_model import EMModel, EMModelTrainer
from joint_snv_mix.classification.em.em_posterior import EMPosterior
from joint_snv_mix.classification.em.joint_models.joint_latent_variables import JointLatentVariables
from joint_snv_mix.classification.utils.beta_binomial_map_estimators import get_mle_p
from joint_snv_mix.classification.utils.log_pdf import log_gamma_pdf, log_beta_pdf, log_beta_binomial_likelihood
import multiprocessing
import numpy as np



class JointBetaBinomialModel( EMModel ):
    def __init__( self ):
        self.trainer_class = JointBetaBinomialModelTrainer
        
        self.log_likelihood_func = joint_beta_binomial_log_likelihood

class JointBetaBinomialModelTrainer( EMModelTrainer ):
    def _init_components( self ):
        self.latent_variables = JointBetaBinomialLatentVariables( self.data )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = JointBetaBinomialPosterior( self.data, self.priors, self.responsibilities )
        
        self.lower_bound = JointBetaBinomialLowerBound( self.data, self.priors )

class JointBetaBinomialLatentVariables( JointLatentVariables ):
    def __init__( self, data ):
        JointLatentVariables.__init__( self, data )
        
        self.likelihood_func = joint_beta_binomial_log_likelihood

class JointBetaBinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_beta_binomial_log_likelihood
        
        self.ncomponents = 9
    
    def _get_log_density_parameters_prior( self ):
        precision_term = 0
        location_term = 0
        
        for genome in range( 2 ):
            for component in range( self.ncomponents ):
                alpha = self.parameters['alpha'][genome][component]
                beta = self.parameters['beta'][genome][component]
                
                s = alpha + beta
                mu = alpha / s
        
                precision_priors = self.priors['precision'][genome][component]
                location_priors = self.priors['location'][genome][component]
                
                precision_term += log_gamma_pdf( s, precision_priors[0], precision_priors[1] )
                location_term += log_beta_pdf( mu, location_priors[0], location_priors[1] )
                
        return precision_term + location_term
        
class JointBetaBinomialPosterior( EMPosterior ):
    def _init_parameters( self ):
        '''
        Initialise parameters using method of moments (MOM) estiamtes.
        '''
        self.ncomponents = 9
        self.p = multiprocessing.Pool( maxtasksperchild=1 )
        
        self.parameters = {}
        
        self._update_mix_weights()
        
#        s = self.priors['precision'][:, :, 0] * self.priors['precision'][:, :, 1] 
#        mu = self.priors['location'][:, :, 0] / \
#            ( self.priors['location'][:, :, 0] + self.priors['location'][:, :, 1] ) 
#        
#        self.parameters['alpha'] = s * mu
#        self.parameters['beta'] = s * ( 1 - mu )

        self.parameters['alpha'] = np.array( [
                                             [99, 5, 1, 99, 5, 1, 99, 5, 1],
                                             [99, 99, 99, 5, 5, 5, 1, 1, 1],
                                             ], dtype=np.float64 )
        self.parameters['beta'] = np.array( [
                                             [1, 5, 99, 1, 5, 99, 1, 5, 99],
                                             [1, 1, 1, 5, 5, 5, 99, 99, 99],
                                             ], dtype=np.float64 )
        
        print "Initial parameter values : ", self.parameters
    
    def _update_density_parameters( self ):        
        print "Begining numerical optimisation of alpha and beta."
        
        vars = []
        
        for genome in range( 2 ):
            for component in range( self.ncomponents ):                
                a = self.data.a[genome]
                b = self.data.b[genome]
                
                x = np.zeros( ( 2, ) )
                
                x[0] = self.parameters['alpha'][genome][component]
                x[1] = self.parameters['beta'][genome][component]
                
                resp = self.responsibilities[:, component]
                
                precision_prior = self.priors['precision'][genome][component]
                location_prior = self.priors['location'][genome][component]
                
                vars.append( [x, a, b, resp, component, precision_prior, location_prior] )
                       
        results = self.p.map( get_mle_p, vars )
        
        for genome in range( 2 ):
            for component in range( self.ncomponents ):
                i = genome * self.ncomponents + component
                
                self.parameters['alpha'][genome][component] = results[i][0]
                self.parameters['beta'][genome][component] = results[i][1]
        
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

    log_likelihoods = normal_log_likelihoods + tumour_log_likelihoods

    pi = parameters['pi']
    log_likelihoods = log_likelihoods + np.log( pi )
    
    return log_likelihoods
