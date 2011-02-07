import numpy as np

from scipy.special import gammaln, betaln

import matplotlib.pyplot as plot

from joint_snv_mix.file_formats.jsm import JointSnvMixReader
from joint_snv_mix import constants

def main( jsm_file_name, prefix ):
    reader = JointSnvMixReader( jsm_file_name )    
    parameters = reader.get_parameters()    
    reader.close()
    
    n = 100
    
    global a, b 
    a = np.linspace( 0, n, n + 1 )
    b = n - a
    
    a = a.reshape( ( a.size, 1 ) )
    b = b.reshape( ( b.size, 1 ) )

    recursive_plot( parameters, prefix )

def recursive_plot( params, plot_name ):    
    for name, value in sorted( params.iteritems() ):
        if name in constants.genomes:
            plot_name = plot_name + "_" + name + ".pdf"
            
            resp = get_marginal_resp( params['pi'], name )            
            
            plot_parameters( value['alpha'], value['beta'], resp, plot_name )
        
        elif isinstance( value, dict ):
            plot_name = plot_name + "_" + name
            
            recursive_plot( value, plot_name )            
        else:
            continue
                
def get_marginal_resp( pi, genome ):
    pi = pi.reshape( ( 3, 3 ) )
    
    if genome == "normal":    
        resp = pi.sum( axis=1 ).reshape( ( 1, 3 ) )
    else:
        resp = pi.sum( axis=1 ).reshape( ( 1, 3 ) )
        
    return resp

def plot_parameters( alpha, beta, resp, plot_name ):
    alpha = alpha.reshape( ( 1, 3 ) )
    beta = beta.reshape( ( 1, 3 ) )
    
    nclass = resp.size
        
    posterior = bb_pdf( a, b, alpha, beta )
    posterior = posterior * resp
    posterior = posterior / posterior.sum( axis=1 ).reshape( ( a.size, 1 ) )

    fig = plot.figure()
    ax = plot.subplot( 111 )
    
    lines = []
    
    for i in range( 3 ):
        p = ax.plot( a, posterior[:, i] )
        lines.append( p )    
            
    line_labels = ['AA', 'AB', 'BB']
    
    fig.legend( lines, line_labels, loc='center right' )    
    
    plot.savefig( plot_name )    

def log_factorial( x ):
    return gammaln( x + 1 )

def log_n_choose_k( k, n ):
    return log_factorial( n ) - log_factorial( n - k ) - log_factorial( k )
    
def bb_pdf( a, b, alpha, beta ):
    n = a + b
    log_pdf = log_n_choose_k( a, n ) + betaln( a + alpha, b + beta ) - betaln( alpha, beta )
    
    return np.exp( log_pdf )
    
def joint_bb_pdf( a_1, a_2, b_1, b_2, alpha_1, alpha_2, beta_1, beta_2 ):
    pdf_1 = bb_pdf( a_1, b_1, alpha_1, beta_1 )
    pdf_2 = bb_pdf( a_2, b_2, alpha_2, beta_2 )
    
    return np.exp( np.log( pdf_1 ) + np.log( pdf_2 ) )

if __name__ == "__main__":
    import sys
    
    jsm_file_name = sys.argv[1]
    prefix = sys.argv[1]
    
    main( jsm_file_name, prefix )
