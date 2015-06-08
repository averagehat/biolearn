from __future__ import division
from utils import min_matrix_ij, doublematrix, emptygraph, \
    parseM, drawgraph, pvar
import numpy as np
from utils import adj_str, test_graph, write_graph
import networkx as nx


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
    Q = np.ma.masked_array(np.zeros(D.shape), mask=False)
    return NJ(D, G, N, Q)


#lost to np.inf -> became -infinity confused me
def _Q(D, Q, n):
    ''' (n-2) * Di - sum(di) - sum(dj)
    transform a distance matrix into neighbor-joining matrix
    which allows minimum to reflect maximum distance from other clusters.'''
    #Q = np.full(D.shape, np.inf)
    #(D[i] * (n - 2)) - ( D[i].sum()  + D.sum(axis=0) )
    #Can't change D here, we need it later!
    shp = D.shape[0]
    '''D.sum(axis=0) sums "down the rows"'''
    ''' add the row-wise sum three times to the column-wise sum.'''
    intermediate = np.tile(D.sum(axis=1), shp).reshape(shp, shp).T + D.sum(axis=0)
    Q[:] = D*(n-2) - intermediate
    Q.mask[:] = D.mask
    return Q
#    def _transform(i):
#         #return D[i] + (D[i] * (n - 2)) - ( D[i].sum()  + D.sum(axis=0) )
#         return D[i]*(n-2) - (D[i].sum() + D.sum(axis=0))
#    Q[:] = map(_transform, xrange(len(Q)))

'''
Neighbor-joining works by iterative clustering. Steadily reduce the size of the matrix by
resolving parent nodes. The tree is constructed from the bottom up, with each iteration (excluding
the base case) forming a new subtree of one parent and two child nodes. The child nodes are then
removed from the matrix and replaced with their parent, which represents the cluster they were absorbed into.
This process is repeated until the number of clusters is 2. At this point the algorithm returns the subtree of
size two with the limb-length representing the actual distance between the clusters. One of these nodes will become the root.'''
'''Instead of reshaping the array or broadcasting a view, we mask. This allows us to avoid expensive concatenation operation
at the expense of creating a larger matrix to begin with. It also simplifies the process of keeping track of clusters; We can
safely refer to each cluster by its index in the matrix, because new clusters are always added to the matrix n+i, where n is the number of
starting leaves and i is the number of iterations. '''
def NJ(D, G, n, Q=None):
    ''' base case '''
    if n == 2:
        #just get non-infinite, argmin() skips masked vals
        #i and j here are the only remaining elements of the matrix.
        i, j = min_matrix_ij(D)
        G.add_edge(i, j, {'weight' : D[i,j]})
        return G
    '''Create our 'enhanced' distance matrix, which gives us information on the overall
    distance of the minimum to other clusters.'''
    Q = _Q(D, Q, n)
    ''' the closest new clusters are the minimum elements of the matrix. '''
    i, j = min_matrix_ij(Q)
    m = len(G)
    add_to_graph(G, D, m, i, j, n)
    ''' the order of the computations and unmasking is important.
    We cannot unmask D[m] before computing u.'''
    u = ( D[i] + D[j] - D[i,j] ) / 2.0
    D.mask[m] = D.mask[:, m] = False
    D[m] = D[:, m] =  u
    ''' remove clusters i and j from the array. '''
    D.mask[i] = D.mask[j] = D.mask[:, i] = D.mask[:, j] = True
    return NJ(D, G, n-1, Q)

wikim=parseM('''
0   5   9   9   8
5   0   10  10  9
9   10  0   8   7
9   10  8   0   3
8   9   7   3   0''')
_G = neighborjoin(wikim)

#drawgraph(_G)

rlm = parseM('''
0   295 300 524 1077    1080    978 941 940
295 0   314 487 1071    1088    1010    963 966
300 314 0   472 1085    1088    1025    965 956
524 487 472 0   1101    1099    1021    962 965
1076    1070    1085    1101    0   818 1053    1057    1054
1082    1088    1088    1098    818 0   1070    1085    1080
976 1011    1025    1021    1053    1070    0   963 961
941 963 965 962 1057    1085    963 0   16
940 966 956 965 1054    1080    961 16  0'''   )
G = neighborjoin(rlm)
animals = 'Cow Pig Horse   Mouse   Dog Cat Turkey  Civet   Human'.split()
mapping =  dict(enumerate(animals))
G = nx.relabel_nodes(G,     mapping)
drawgraph(G, edgekey=False, big=True)


exp = '''
0->4:8.000
1->5:13.500
2->5:16.500
3->4:12.000
4->5:2.000
4->0:8.000
4->3:12.000
5->1:13.500
5->2:16.500
5->4:2.000
'''.strip()

excer=parseM('''0   23  27  20
23  0   30  28
27  30  0   30
20  28  30  0''')
AG = neighborjoin(excer)
actual = adj_str(AG, as_float=True).strip()
assert sorted(exp) == sorted(actual)
print "passed!"

test_graph('datas/NJ_input.txt', 'datas/nj_exp.txt', neighborjoin)
#write_graph('datas/njin.txt', 'datas/nj.txt', neighborjoin)
