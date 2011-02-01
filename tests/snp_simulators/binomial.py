'''
Created on 2010-11-18

@author: Andrew Roth
'''
import numpy as np

from scipy.stats import beta, binom, poisson

from tables import *

def draw_joint_sample( params, n, mean_depth ):    
    class_params = []
    
    for i in range( 3 ):
        for j in range( 3 ):
            class_params.append( [ 
                                 params['mu_1'][i], params['mu_2'][j] 
                                 ] )
    
    multinomial_draw = np.random.multinomial( 1, params['pi'], n )

    labels = np.argmax( multinomial_draw == 1, axis=1 )
    
    a_b = []

    for label in labels:
        a_b.append( class_params[label] )
    
    a_b = np.array( a_b )

    x_N = draw_binomial_sample( a_b[:, 0], mean_depth, n )
    x_T = draw_binomial_sample( a_b[:, 1], mean_depth, n )
        
    x = np.hstack( ( x_N, x_T ) )
    
    return x, labels

def draw_binomial_sample( mu, mean_depth, n ):
    depth_draw = poisson.rvs( mean_depth, size=n )
    
    ref_match = binom.rvs( depth_draw, mu )
    
    x = np.column_stack( ( ref_match, depth_draw ) )
    
    return x

def write_counts( file_name, counts, labels, params ):
    filter = Filters( complevel=1, complib='zlib' )
    
    h5_file = openFile( file_name, 'w', filters=filter )
    
    atom = UInt16Atom( () )
    shape = counts.shape
    
    h5_file.createCArray( h5_file.root, 'counts', atom, shape )
    h5_file.root.counts[:] = counts[:]
    
    atom = UInt8Atom( () )
    shape = labels.shape
    
    h5_file.createCArray( h5_file.root, 'labels', atom, shape )
    h5_file.root.labels[:] = labels[:]
    
    params_group = h5_file.createGroup( h5_file.root, 'parameters' )
    
    shape = params['mu_1'].shape
    atom = Float64Atom( () )
    
    h5_file.createCArray( params_group, 'mu_1', atom, shape )
    h5_file.createCArray( params_group, 'mu_2', atom, shape )
    
    params_group.mu_1[:] = params['mu_1'][:]
    params_group.mu_2[:] = params['mu_2'][:]
    
    shape = params['pi'].shape
    
    h5_file.createCArray( params_group, 'pi', atom, shape )
    
    params_group.pi[:] = params['pi'][:]
    
    h5_file.close()

if __name__ == "__main__":
    n = int( 1e4 )
    mean_depth = 20
    
    params = {}
    
    params['pi'] = np.array( [1e3, 10, 5, 5, 100, 5, 1, 2, 100] )
    params['pi'] = params['pi'] / params['pi'].sum()
    
    
    params['mu_1'] = np.array( [0.99, 0.5, 0.01] )
    params['mu_2'] = np.array( [0.99, 0.5, 0.01] )
    
    counts, labels = draw_joint_sample( params, n, mean_depth )
    
    write_counts( '../data/binomial_small.sim', counts, labels, params )
