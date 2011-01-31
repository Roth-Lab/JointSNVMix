import numpy as np

from scipy.special import gammaln, betaln

import matplotlib.pyplot as plot

from joint_snv_mix.file_formats.jsm import JointSnvMixReader

def main( jsm_file_name ):
    reader = JointSnvMixReader( jsm_file_name )
    
    parameters = reader.get_parameters()
    
    reader.close()
    
    n = 100
    
    a = np.linspace( 0, n, n + 1 )
    b = n - a
    
    alpha = parameters['alpha'][0].reshape( ( 3, 1 ) )
    beta = parameters['beta'][0].reshape( ( 3, 1 ) )
    
    pi = parameters['pi'].reshape( ( 3, 3 ) )
    
    mw = pi.sum( axis=1 ).reshape( ( 3,1 ) )
       
    p = bb_pdf( a, b, alpha, beta )
    p = p * mw
    p = p / p.sum( axis=0 ).reshape( ( 1, p.shape[1] ) )

    plot.plot( a, p[0] )
    plot.plot( a, p[1] )
    plot.plot( a, p[2] )
    
    plot.savefig( 'dist.pdf' )

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
    
    main( jsm_file_name )
