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
                                 params['alpha'][0, i], params['beta'][0, i],
                                 params['alpha'][1, j], params['beta'][1, j] 
                                 ] )
    multinomial_draw = np.random.multinomial( 1, params['pi'], n )

    labels = np.argmax( multinomial_draw == 1, axis=1 )
    
    a_b = []

    for label in labels:
        a_b.append( class_params[label] )
    
    a_b = np.array( a_b )

    x_N = draw_bb_sample( a_b[:, 0], a_b[:, 1], mean_depth, n )
    x_T = draw_bb_sample( a_b[:, 2], a_b[:, 3], mean_depth, n )
        
    x = np.hstack( ( x_N, x_T ) )
    
    return x, labels

def draw_bb_sample( a, b, mean_depth, n ):
    mu = beta.rvs( a, b )
    
    # Do this to avoid zero depth and still keep mean depth.
    d = poisson.rvs( mean_depth - 1, size=n ) + 1
    
    a = binom.rvs( d, mu )
    
    b = d - a
    
    x = np.column_stack( ( a, b ) )
    
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
    
    shape = params['alpha'].shape
    atom = Float64Atom( () )
    
    h5_file.createCArray( params_group, 'alpha', atom, shape )
    h5_file.createCArray( params_group, 'beta', atom, shape )
    
    params_group.alpha[:] = params['alpha'][:]
    params_group.beta[:] = params['beta'][:]
    
    shape = params['pi'].shape
    
    h5_file.createCArray( params_group, 'pi', atom, shape )
    
    params_group.pi[:] = params['pi'][:]
    
    h5_file.close()

if __name__ == "__main__":
    n = int( 1e5 )
    mean_depth = 40
    
    params = {}
    
    params['pi'] = np.array( [1e4, 10, 5, 4, 100, 4, 0.1, 0.2, 100] )
    params['pi'] = params['pi'] / params['pi'].sum()
    
    
    params['alpha'] = np.array( [
                                 [999., 2., 5],
                                 [999., 10., 0.9]
                                ] )
    
    params['beta'] = np.array( [
                                 [10, 2., 80.],
                                 [0.99, 10., 120.] 
                                 ] )
    
    counts, labels = draw_joint_sample( params, n, mean_depth )
    
    write_counts( '../data/bb_40.sim', counts, labels, params )
