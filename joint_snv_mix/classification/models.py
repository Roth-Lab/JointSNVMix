'''
Created on 2010-11-22

@author: Andrew Roth
'''
import numpy as np
#np.seterr( invalid='raise' )

from joint_snv_mix.classification.latent_variables import IndependentBinomialLatentVariables, IndependentBetaBinomialLatentVariables, JointBetaBinomialLatentVariables, JointBinomialLatentVariables
from joint_snv_mix.classification.likelihoods import independent_binomial_log_likelihood, independent_beta_binomial_log_likelihood, joint_beta_binomial_log_likelihood, joint_binomial_log_likelihood
from joint_snv_mix.classification.lower_bounds import IndependenBinomialLowerBound, IndependentBetaBinomialLowerBound, JointBetaBinomialLowerBound, JointBinomialLowerBound
from joint_snv_mix.classification.posteriors import IndependentBinomialPosterior, IndependentBetaBinomialPosterior, JointBetaBinomialPosterior, JointBinomialPosterior
from joint_snv_mix.classification.utils.normalise import log_space_normalise_rows

class EMModel( object ):
    def __init__( self ):
        self.trainer_class = None
        self.log_likelihood_func = None
        
    
    def train( self, data, priors, max_iters, tolerance ):
        '''
        Train the model using EM.
        
        Input: JointData object
        '''   
        trainer = self.trainer_class( data, max_iters, tolerance, priors )
        
        parameters = trainer.run()
                
        return parameters

    def classify( self, data, parameters ):
        log_responsibilities = self.log_likelihood_func( data, parameters )
        
        responsibilities = log_space_normalise_rows( log_responsibilities )
        
        return responsibilities

class EMModelTrainer( object ):
    def __init__( self, data, max_iters, tolerance, priors ):
        self.max_iters = max_iters
        
        self.tolerance = tolerance
        
        self.data = data

        self.priors = priors
            
        self._init_components()
        
    def run( self ):
        iters = 0
        converged = False
        
        parameters = self.posterior.parameters
        old_posterior_value = self.lower_bound.get_lower_bound( parameters )
  
        while not converged:
            self._M_step()
            self._E_step()

            posterior_value = self.lower_bound.get_lower_bound( self.parameters )
            
            if iters > 0:
                posterior_change = ( posterior_value - old_posterior_value ) / abs( old_posterior_value )
            else:
                posterior_change = float( 'inf' )
            
            self._print_diagnostic_message( iters, posterior_value, old_posterior_value, posterior_change )
            old_posterior_value = posterior_value
            
            if posterior_change < 0:
                print "Posterior decreased. This could be a bug or overly stringent convergence criterion."
                converged = True
            
            elif posterior_change < self.tolerance:
                converged = True
                            
            if iters >= self.max_iters:
                print "Maximum numbers of EM iterations exceeded. Exiting training."                
                converged = True
            
            iters += 1
                     
        return self.parameters
                  
    def _E_step( self ):
        self.latent_variables.update( self.parameters )
        self.responsibilities = self.latent_variables.responsibilities
        
    def _M_step( self ):
        self.posterior.update( self.responsibilities )
        self.parameters = self.posterior.parameters
        
    def _print_diagnostic_message( self, iters, posterior_value, old_posterior_value, posterior_change ):
        print "#" * 100
        print "# Diagnostics."
        print "#" * 100
        print "Number of iterations : ", iters
        print "New posterior : ", posterior_value
        print "Old posterior : ", old_posterior_value 
        print "Posterior change : ", posterior_change
    
        print "Parameters :"
        
        for param_name, param_value in self.posterior.parameters.items():
            print param_name, param_value     
    
    def _init_components( self ):
        raise NotImplemented

#=======================================================================================================================
# Independent Models
#=======================================================================================================================
class IndependenBetaBinomialModel( EMModel ):
    def __init__( self ):
        self.trainer_class = IndependenBetaBinomialTrainer
        
        self.log_likelihood_func = independent_beta_binomial_log_likelihood

class IndependenBetaBinomialTrainer( EMModelTrainer ):
    def _init_components( self ):
        self.latent_variables = IndependentBetaBinomialLatentVariables( self.data )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = IndependentBetaBinomialPosterior( self.data, self.priors, self.responsibilities )
        
        self.lower_bound = IndependentBetaBinomialLowerBound( self.data, self.priors )
        
class IndependentBinomialModel( EMModel ):
    def __init__( self ):
        self.trainer_class = IndependentBinomialModelTrainer
        
        self.log_likelihood_func = independent_binomial_log_likelihood

class IndependentBinomialModelTrainer( EMModelTrainer ):
    def _init_components( self ):
        self.latent_variables = IndependentBinomialLatentVariables( self.data )
        
        self.responsibilities = self.latent_variables.responsibilities
        
        self.posterior = IndependentBinomialPosterior( self.data, self.priors, self.responsibilities )
        
        self.lower_bound = IndependenBinomialLowerBound( self.data, self.priors )
        
#=======================================================================================================================
# Joint Models
#=======================================================================================================================
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
