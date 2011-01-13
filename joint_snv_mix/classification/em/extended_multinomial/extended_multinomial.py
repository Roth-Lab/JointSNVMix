'''
Created on 2010-12-09

@author: Andrew Roth
'''
import numpy as np

from joint_snv_mix.classification.em.em_lower_bound import EMLowerBound
from joint_snv_mix.classification.em.em_model import EMModel, EMModelTrainer
from joint_snv_mix.classification.em.em_posterior import EMPosterior
from joint_snv_mix.classification.em.extended_multinomial.extended_multinomial_latent_variables import ExtendedMultinomialLatentVariables
from joint_snv_mix.classification.utils.log_pdf import log_dirichlet_pdf, log_multinomial_likelihood

nclass = 15

class JointExtendedMultinomialModel( EMModel ):
    def __init__( self ):
        self.trainer_class = JointExtendedMultinomialModelTrainer
        
        self.log_likelihood_func = joint_multinomial_log_likelihood

class JointExtendedMultinomialModelTrainer( EMModelTrainer ):
    def _init_components( self ):
        self.latent_variables = JointExtendedMultinomialLatentVariables( self.data )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = JointExtendedMultinomialPosterior( self.data, self.priors, self.responsibilities )
        
        self.lower_bound = JointExtendedMultinomialLowerBound( self.data, self.priors )

class JointExtendedMultinomialLowerBound( EMLowerBound ):
    def __init__( self, data, priors ):  
        EMLowerBound.__init__( self, data, priors )
        
        self.log_likelihood_func = joint_multinomial_log_likelihood
    
    def _get_log_density_parameters_prior( self ):
        log_prior = 0.
        
        for genome in range( 2 ):
            rho = self.parameters['rho'][genome]
            delta = self.priors['delta'][genome]
            
            log_prior += np.sum( log_dirichlet_pdf( rho, delta ) )
        
        log_prior = log_prior.sum()

        return log_prior
    
class JointExtendedMultinomialLatentVariables( ExtendedMultinomialLatentVariables ):
    def __init__( self, data ):
        ExtendedMultinomialLatentVariables.__init__( self, data )
        
        self.likelihood_func = joint_multinomial_log_likelihood
        
class JointExtendedMultinomialPosterior( EMPosterior ):
    def _init_parameters( self ):
        self.parameters = {}
        
        self._update_mix_weights()
        
        self.parameters['rho'] = []
        
        for i in range( 2 ):
            delta = self.priors['delta'][i]
            
            shape = ( delta.shape[0], 1 )
            
            rho = delta / delta.sum( axis=1 ).reshape( shape )
            
            self.parameters['rho'].append( rho )
           
    def _update_density_parameters( self ):
        marginals = get_marginals( self.responsibilities )
        
        for genome in range( 2 ):
            counts = self.data.counts[genome]
            
            delta = self.priors['delta'][genome]
            
            tau = marginals[genome]
            
            self.parameters['rho'][genome] = self._update_rho( counts, tau, delta )
    
    def _update_rho( self, counts, tau, delta ):       
        counts = counts.reshape( ( 1, counts.shape[0], counts.shape[1] ) )
        
        tau = np.swapaxes( tau, 0, 1 )
        
        tau = tau.reshape( ( tau.shape[0], tau.shape[1], 1 ) )
        
        marginal_counts = tau * counts
        
        marginal_counts = marginal_counts.sum( axis=1 )
        
        numerator = marginal_counts + delta - 1
        
        denominator = numerator.sum( axis=1 ).reshape( ( numerator.shape[0], 1 ) )
        
        return np.exp( np.log( numerator ) - np.log( denominator ) )
        
def joint_multinomial_log_likelihood( data, parameters ):
    counts_1 = data.counts[0]
    counts_2 = data.counts[1]
    
    rho_1 = parameters['rho'][0]
    rho_2 = parameters['rho'][1]

    
    normal_log_likelihoods = log_multinomial_likelihood( counts_1, rho_1 )
    tumour_log_likelihoods = log_multinomial_likelihood( counts_2, rho_2 )

    column_shape = ( normal_log_likelihoods[:, 0].size, 1 )

    log_likelihoods = []
    
    for i in range( nclass ):
        log_likelihoods.append( normal_log_likelihoods[:, i].reshape( column_shape ) + tumour_log_likelihoods )

    log_likelihoods = np.hstack( log_likelihoods )

    pi = parameters['pi']
    log_pi = np.log( pi )

    log_likelihoods = log_likelihoods + log_pi
    
    return log_likelihoods

def get_marginals( responsibilities ):
    marginals = []

    marginals.append( get_normal_marginals( responsibilities ) )
    marginals.append( get_tumour_marginals( responsibilities ) )

    return marginals

def get_normal_marginals( responsibilities ):
    nrows = responsibilities.shape[0]

    shape = ( nrows, nclass )

    normal_marginals = np.zeros( shape )
    
    for i in range( nclass ):
        start = nclass * i
        stop = nclass * i + nclass
        indices = range( start, stop )
        normal_marginals[:, i] = np.sum( responsibilities[:, indices], axis=1 )

    return normal_marginals

def get_tumour_marginals( responsibilities ):
    nrows = responsibilities.shape[0]

    shape = ( nrows, nclass )

    tumour_marginals = np.zeros( shape )
    
    for i in range( nclass ):
        indices = range( i, nclass ** 2, nclass )

        tumour_marginals[:, i] = np.sum( responsibilities[:, indices], axis=1 )
    
    return tumour_marginals
