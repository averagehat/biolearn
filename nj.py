from __future__ import division
from utils import min_matrix_ij, doublematrix, emptygraph, \
    parseM, drawgraph, pvar
import numpy as np

def get_limb_weights(D, i, j, n):
    W_i = (D[i,j])/2.0  +  (1/(2.0*(n-2)) * (D[i].sum() - D[j].sum()))
    W_j = D[i, j] - W_i
    return W_i, W_j

def add_to_graph(G, D, m, c1, c2, n):
    ''' m is the new parent '''
    C1_w, C2_w = get_limb_weights(D, c1, c2, n)
    G.add_node(m)
    G.add_edge(m, c1, {'weight' : C1_w})
    G.add_edge(m, c2, {'weight' : C2_w})
    return G

def neighborjoin(M):
    N = M.shape[0]
    D = np.zeros((N*2, N*2))
    #D = np.zeros((N*2, N*2), np.nan)
    D[0:N, 0:N] = M
    D = np.ma.masked_array(D, mask=False)
    D.mask[N:] = D.mask[:, N:] = True
    D.mask[np.diag_indices(len(D))] = True
    G = emptygraph(N)
    return NJ(D, G, N)


#lost to np.inf -> became -infinity confused me
def _Q(D, n):
    ''' (n-2) * Di - sum(di) - sum(dj)
    transform a distance matrix into neighbor-joining matrix
    which allows minimum to reflect maximum distance from other clusters.'''
    #Q = np.full(D.shape, np.inf)
    Q = D.copy()

    #(D[i] * (n - 2)) - ( D[i].sum()  + D.sum(axis=0) )
    #Ack! can't change D here, we need it later!
    def _transform(i):
         #return D[i] + (D[i] * (n - 2)) - ( D[i].sum()  + D.sum(axis=0) )
         return D[0]*(n-2) - (D[0].sum() + D.sum(axis=0))
    Q[:] = map(_transform, xrange(len(Q)))
    Q.mask[:] = D.mask

    #Q = np.apply_along_axis(_transform, 0, D)
#    for i, _ in enumerate(D):
#        Q[i] = (D[i] * (n - 2)) - ( D[i].sum()  + D.sum(axis=0) )
#        D[i] *= (n - 2)
#        #D[i] -= D[i].sum()
#        #D[i] -= D.sum(axis=0)
#        D[i] -= (D[i].sum() + D.sum(axis=0))
    return Q

def NJ(D, G, n=None):
    n = len(D)//2 if n is None else n
    if n == 2:
        #just get non-infinite, argmin() skips masked vals
        i, j = min_matrix_ij(D)
        G.add_edge(i, j, {'weight' : D[i,j]})
        return G
    Q = _Q(D, n)
    print D
    print Q
    i, j = min_matrix_ij(Q)
    pvar(i), pvar(j)
    m = len(G)
    D[m] = D[:, m] = ( D[i] + D[j] - D[i,j] ) / 2.0
    add_to_graph(G, D, m, i, j, n)
    #D.mask[m] = False
    #Dm = phyloavg
    D.mask[i] = D.mask[j] = D.mask[:, i] = D.mask[:, j] = True
    return NJ(D, G, n-1)

wikim=parseM('''
0   5   9   9   8
5   0   10  10  9
9   10  0   8   7
9   10  8   0   3
8   9   7   3   0''')

G = neighborjoin(wikim)

drawgraph(G)

