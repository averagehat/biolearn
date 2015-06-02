import numpy as np
from fn import F
import networkx as nx
from itertools import repeat
from scrach import parseM

def min_matrix_ij(M):
    ''' :return indices (i, j) for the minimum element in the matrix.
    assumes diagonals != 0 or are masked. '''
    res= np.unravel_index(M.argmin(), M.shape)
    return res

''' the closest cluster is the minimum elemnt in the *distance* matrix!'''
closest_clusters = min_matrix_ij

'''
Learned: read the whole algorithm (or stick to one?)
original implementation had ages wrong because
#NOTE: recursion is nice for debugging up the stack-frame. (if pass immutable references, like M.copy())
'''
def UPGMA(m):
    '''number of clusters = N (in an NxN matrix).'''
    '''First, create a double-sized matrix to make it easy to replace clusters.
    Actually, should just right over one of the rows at random? I think. '''
    M = np.full((n*2, n*2), np.inf)
    M[0:n, 0:n] = m
    ''' the mask is needed to skip diagonals.'''
    #TODO: get indices from min_matrix without masking!
    #M = np.ma.masked_array(M, mask=False)
    #M.mask[np.diag_indices(len(M))] = True
    ''' nvm, this skips diagonals, because infinity will never be the minimum.'''
    M[np.diag_indices(len(M))] = np.inf
    N = M.shape[0]/2
    G = nx.DiGraph()
    ''' add blank nodes to ease the construction of the graph.'''
    G.add_nodes_from(zip(xrange(N), [dict(weight=0) for i in xrange(N)]))
    agemap={}
    for i in xrange(N):
        agemap[i] = 0
    G2, agemap= _UPGMA(M.copy(), G, N, {}, agemap)
    ''' set the weights so that the limbs are not too long.
    Note that the average distance (ll) computed between clusters
    is the universal distance (from the original (now clustered) leaves!'''
    for k, v in G2.edges():
        G2.get_edge_data(k, v)['weight'] = agemap[k] - agemap[v]
    return G2

def clustdist(D, i, j, wi, wj, c):
    '''i, j are clusters, (wi, wj) = weight of i and j'''
    return (D[i, c]*wi + D[j, c]*wj)/(wi + wj)

def allbut(end, excluding):
    ''' get range(0, end) excluding one element '''
    if excluding < end: return range(start or 0, excluding) + range(excluding + 1, end)
    else: return range(0, end)

def _UPGMA(M, G, n, clustweight, agemap):
    ''' M should be size n*2 and half-masked (or infinity? infinity never gets picked).'''
    if n == 1: return   G, agemap
    c1, c2 = closest_clusters(M)
    #typos: big error here where c1 was repeated. ie w1, w2 = get(c1), get(c1)!
    w1, w2 = clustweight.get(c1, 1), clustweight.get(c2, 1)

    '''a cluster is defined by its mean, so their root seperates leaves evenly.'''
    ll = M[c1, c2]/2.0
    #print M #mildly useful
    '''give C a random name but small enough to place within the matrix. works because max nodes < M.shape'''
    C = len(G)
    ''' clusters are weight by the amount of leaves they represent!'''
    clustweight[C] = w1 + w2
    G.add_node(C)
    G.add_edge(C, c1, dict(weight=0))
    G.add_edge(C, c2, dict(weight=0))
    '''save for later'''
    agemap[C] = ll
    othercs = xrange(len(M))
    _dist = F(clustdist, M, c1, c2, w1, w2)
    #TODO: better way to do this, because don't actually have to ignore diagonals
    M[C] = M[:, C] = map(_dist, othercs)
    ''' remove derived clusters c1 and c2 by setting them to infinity; this works
    because our algorithm seeks minimum i think'''
    M[c1] = M[c2] = M[:, c1] = M[:, c2] = np.inf
    return _UPGMA(M.copy(), G, n-1, clustweight, agemap)

m=parseM('''
0   20  17  11
20  0   20  13
17  20  0   10
11  13  10  0''')
n = len(m)
G = UPGMA(m)
    #M.mask[c1] = M.mask[c2] = M.mask[:, c1] = M.mask[:, c2] = np.inf
    #M.mask[C] = M.mask[:, C] = False
    #M.mask[c1] = M.mask[c2] = M.mask[:, c1] = M.mask[:, c2] = True
