'''
Created on 2010-10-14

@author: Andrew Roth
'''
import numpy as np

def get_marginals( responsibilities ):
    marginals = []

    marginals.append( get_normal_marginals( responsibilities ) )
    marginals.append( get_tumour_marginals( responsibilities ) )

    return marginals

def get_normal_marginals( responsibilities ):
    nrows = responsibilities.shape[0]

    shape = ( nrows, 3 )

    normal_marginals = np.zeros( shape )

    normal_marginals[:, 0] = responsibilities[:, 0] + responsibilities[:, 1] + responsibilities[:, 2]
    normal_marginals[:, 1] = responsibilities[:, 3] + responsibilities[:, 4] + responsibilities[:, 5]
    normal_marginals[:, 2] = responsibilities[:, 6] + responsibilities[:, 7] + responsibilities[:, 8]

    return normal_marginals

def get_tumour_marginals( responsibilities ):
    nrows = responsibilities.shape[0]

    shape = ( nrows, 3 )

    tumour_marginals = np.zeros( shape )

    tumour_marginals[:, 0] = responsibilities[:, 0] + responsibilities[:, 3] + responsibilities[:, 6]
    tumour_marginals[:, 1] = responsibilities[:, 1] + responsibilities[:, 4] + responsibilities[:, 7]
    tumour_marginals[:, 2] = responsibilities[:, 2] + responsibilities[:, 5] + responsibilities[:, 8]

    return tumour_marginals
