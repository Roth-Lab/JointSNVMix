import numpy as np

from scipy.special import gammaln, betaln

import matplotlib.pyplot as plot

from joint_snv_mix.file_formats.jsm import JointSnvMixReader

def main( jsm_file_name, out_file_name ):
    reader = JointSnvMixReader( jsm_file_name )
    parameters = reader.get_parameters()
    reader.close()

    n = 50

    a = np.linspace( 0, n, n + 1 )
    b = n - a

#    alpha_N = np.array( [100, 5, 1], dtype = np.float ).reshape( ( 3, 1 ) )
#    beta_N = np.array( [1, 5, 100], dtype = np.float ).reshape( ( 3, 1 ) )
#
#    alpha_T = alpha_N * 20
#    beta_T = beta_N * 20
#
#    pi = np.array( [1e6, 1e2, 10, 1e2, 1e4, 1e2, 1, 5, 1e4] )
#    pi = pi / pi.sum()

    alpha_N = parameters['alpha'][0].reshape( ( 3, 1 ) )
    beta_N = parameters['beta'][0].reshape( ( 3, 1 ) )

    alpha_T = parameters['alpha'][1].reshape( ( 3, 1 ) )
    beta_T = parameters['beta'][1].reshape( ( 3, 1 ) )

    pi = parameters['pi'].reshape( ( 3, 3 ) )

    mw_N = pi.sum( axis = 1 ).reshape( ( 3, 1 ) )
    mw_T = pi.sum( axis = 0 ).reshape( ( 3, 1 ) )

    params = [
              ( alpha_N, beta_N, mw_N, 'normal' ),
              ( alpha_T, beta_T, mw_T, 'tumour' ),
              ]

    fig = plot.figure( figsize = ( 8, 6 ) )

    for i, params in enumerate( params ):
        alpha = params[0]
        beta = params[1]
        mw = params[2]
        sample_name = params[3]

        p = bb_pdf( a, b, alpha, beta )
        p = p * mw
        nc = p.sum( axis = 0 ).reshape( ( 1, p.shape[1] ) )
        p = p / nc

        ax = fig.add_subplot( 2, 1, i + 1 )

        ax.grid( True )
        ax.set_title( sample_name )

        for i in range( len( p ) ):
            ax.plot( a, p[i] )

    plot.savefig( out_file_name )

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
    out_file_name = sys.argv[2]

    main( jsm_file_name, out_file_name )
