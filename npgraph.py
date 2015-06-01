'''
represent visitors as a boolean array
represent tentative distance as float array.
represent pred. as a node (int) array.  

initialize visitors.
initialize preds as NaN.
initialize tentative distance as inf, excepting start node =0

select a current node.
for each of current node's neighbors:
    1. update neighbor's tentative distance & predecessor if it's less.
mark current node as visited.
restart with current node= node with minimum tentative 
'''

import numpy as np
from itertools import ifilterfalse
from operator import methodcaller as mc
from functools import partial
from fn import _ as X, F
import operator as op
def filterfalse(func, seq): return list(ifilterfalse(func, seq))
def allbut(start, skip, end):
    return range(start, skip) + range(skip+1, end)

def nbr_idx(M, N):
  #return M[allbut(0, N, M.shape[0])].nonzero()
  candidates = allbut(0, N, M.shape[0])
   #TODO: find a beter way to filter by a boolean matrix.
  return filterfalse(lambda x: np.isnan(M[N][x]), candidates)

def update(M, N, v):
    M[N] = v; M[:, N] = v

def iterinfo(G, node, dists, vstd):
    ''' return: (idx, tentativedist, weigt) '''
    #TODO: could simply multiply with visited before calling nbr_indx
    for i in nbr_idx(G, node):
        if not vstd[i]: yield i, dists[i], G[node][i]

#TODO: this actually assumes starat is always zero... see [0] + ....
def init(G, start=0):
    ''' :return (visitors, tent. distance, tent. pred) '''
    N = G.shape[0]
    vstd = np.full(N, False)
    dists = np.array([0]+[np.inf]*(N-1))
    preds = np.full(N, np.nan)
    return vstd, dists, preds

def djs(G, start=0):
    '''return pred and distance arrays.'''
    assert (G[np.diag_indices(G.shape[0])] == 0).all(), 'expected zero diagonals.'
    #assert (G.T == G).all(), 'expected symmetric matrix.'
    vstd, dists, preds = init(G, start)
    return _djs(G, start, preds, dists, vstd) 

#TODO: use a mask instead? 
def _djs(G, node, preds, dists, vstd): 
    '''
    for each of current node's neighbors:
        1. update neighbor's tentative distance & predecessor if it's less.
    mark current node as visited.
    restart with current node=node with minimum tentative distance.
    '''
    for idx, t_dist, weight in iterinfo(G, node, dists, vstd):
        walk_dist = dists[node] + weight
        if walk_dist < t_dist:
            preds[idx], dists[idx] = node, walk_dist
    vstd[node] = True
    #unvisited = vstd.nonzero()
    #TODO: find a beter way to filter by a boolean matrix.
    unvisited = filterfalse(vstd.__getitem__, xrange(G.shape[0]))
    if not unvisited:
        return preds, dists
    nextnode = min(unvisited, key=dists.__getitem__)
    return _djs(G, nextnode, preds, dists, vstd)

def unvisited(G, vstd_A):
    return filterfalse(vstd.__getitem__, xrange(G.shape[0]))



'''
if the current node is the target, return current node
else, try path for each unvisited neighobr.
filter on the target-reaching path and return it.
'''
def uniq_path(G,  target, vstd=None, node=0): 
    ''' 
    :param vstd: 1D array
    :return a 2Xd matrix of distances. '''
    if vstd is None:
        vstd = np.full(G.shape[0], False)
    if vstd[node]: return [None]
    if node == target:
        return [node]
    else:
        vstd[node] = True
        nbrs = filterfalse(vstd.__getitem__, nbr_idx(G, node))
        getpath = partial(uniq_path, G, target, vstd)
        path_candidates = filterfalse(mc('__contains__', None), map(getpath, nbrs))
        #path_candidates = filterfalse(F(op.contains, None), map(getpath, nbrs))
        if path_candidates:
            return [node] + path_candidates[0]
        else:
            return [None]






