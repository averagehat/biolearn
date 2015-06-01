import numpy as np
from functools import partial
from operator import itemgetter, ne
from func import compose, compose_all, partial2, ilen, nameddict
from trees import get_products
F = partial
def both(f, g):
    def both_inner(*args, **kwargs):
        first = f(*args, **kwargs)
        return first, g(*args, **kwargs)
    return both_inner

def npGraph(matrix):
    M=np.ma.masked_array(matrix, mask=False, dtype=float)
    get = M.__getitem__
    _set = M.__setitem__
    def update(i, j, value):
        M[i,j]=value
        M[j,i]=value
    def neighborseq(N, val):
        if val == 0:
            raise ValueError
        return (get(N) == val).nonzero()
    def hide_node(N):
        M.mask[N] = M.mask[:, N] = True
    def rowndiag(start): range(0, start) + range(start+1, M.shape[0])
    def colndiag(start):  range(0, start) + range(start+1, M.shape[1])
    def row_non_diag(N): return get((N, rowndiag(N)))
    diag_idx = F(np.diag_indices, M.shape[0])
    diag = compose(get, diag_idx)
    is_neighbor = compose(has_edge, row_non_diag)
    neighbor_weights = compose_all(get, is_neighbor)
    neighbors = compose(np.nonzero, is_neighbor)
    degree = compose(ilen, neighbor_weights)
    update_row = lambda i, v: update(i=i, j=slice(None), val=v)
    update_col = lambda j, v: update(i=slice(None), j=j, val=v)
    #update_row = partial(update, j=slice(None))
    #update_col = partial(update, i=slice(None))
    update_node = both(update_row, update_col)
    update_row_n = partial(update, j=rowndiag)
    update_col_n = partial(update, i=colndiag)
    update_neighbors = both(update_row_n, update_col_n)
    drop_node = partial(update_node, value=NE)
    drop_edge = partial2(_set, NE)
    closest_neighbor = compose_all(M.argmin, np.min, neighbor_weights)
    rowproducts = compose(get_products, rowndiag)
    rowproducts3 = compose(F(get_products, times=3), rowndiag)
    one_neighbor_eq = compose(itemgetter(0), neighborseq)

    return nameddict('NPGraph', {'drop_node' : drop_node, 'drop_edge' : drop_edge, 'closest_neighbor' : closest_neighbor, 'rowproducts3' : rowproducts3,
                      'rowprdocts' : rowproducts, 'one_neighbor_eq' : one_neighbor_eq, 'neighborseq' : neighborseq,
                      'hide_node' : hide_node, 'neighbors' : neighbors, 'degree' : degree, 'neighbor_weights' : neighbor_weights,
                       'update_neighbors' : update_neighbors, 'M' : M, 'is_neighbor' : is_neighbor, 'update_node' : update_node})

NE = np.nan
has_edge = partial(ne, NE)
class Graph(object):
    NE = np.nan

    def __init__(self, M):
        self.M=np.ma.masked_array(M, mask=False, dtype=float)
        self.M[np.diag_indices(M.shape[0])] = 0
        assert np.array_equal(self.M, self.M.T), "Tranpsose matrix not equal"

    def rowndiag(self, start): range(0, start) + range(start+1, self.M.shape[0])
    def colndiag(self, start):  range(0, start) + range(start+1, self.M.shape[1])
    def neighbors(self, N, row=True): return has_edge(self.rowndiag(N)) if row else has_edge(self.colndiag(N))
    def neighbors_idx(self, N, row=True):
        return (self.neighbors(N)).nonzero()
    def neighbor_weights(self, N):
        return self.M[self.neighbors(N)]
    def closest_neighbor(self, N):
        return self.M.argmin(self.neighbor_weights(N).min())
    def neighborseq(self, N, val):
        if val == 0:
            raise ValueError
        return (self.M[N] == val).nonzero()
    def hide_node(self, N):
        self.M.mask[N] = self.M.mask[:, N] = True


def leafpath(self, i, j, neighbors, target, rowwise=True):
    if neighbors[target] != NE
        return neighbors[target]
    else:
        self.hide_node(idx)
    if rowwise:
        next_neighbors = self.neighborseq(row=(not rowwise))
        #return  self.M[i, j] sum([self.leafpath(i, idx, M[:, idx],
'''
create distance matrix

1. empty pred matrix
2. for each node
   1. take path
   2. if current.dist + path_dist < next node's mindist, update it's mindist & predecessor
mark current as visited

pick min visit node (shortest tentative distance), assign as current node

'''







