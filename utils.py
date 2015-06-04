from __future__ import print_function
from func import compose_all, compose, starcompose, dictzip, ilen, pmap
from fn import F, _ as X
from itertools import groupby, starmap, imap, ifilter, izip
from operator import methodcaller as mc, ne
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt

def parseM(raw):
    '''parse & return a space-seperated matrix.'''
    _in = filter(bool, raw.split('\n'))
    return  np.matrix(map(pmap(float), map(str.split, _in)))

def quantify(iterable, pred=bool):
    '''https://docs.python.org/2/library/itertools.html#recipes
    "Count how many times the predicate is true"'''
    return sum(imap(pred, iterable))

def drawgraph(G, edgekey='weight', big=False, **kwargs):
    if big: fig = plt.figure(figsize = (15, 10))
    pos=nx.spring_layout(G)
    nx.draw_networkx(G, pos=pos, **kwargs)
    if edgekey:
        edge_labels=dict([((u,v,),d.get(edgekey, ''))
                         for u,v,d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)#, **kwargs)
    plt.show()
def hamming(s1, s2):
    assert len(s1) == len(s2)
    return ilen(ifilter(bool, imap(ne, s1, s2)))

def logcall(func):
    def wrap(*args, **kwargs):
        print( func.__name__)
        print( args, kwargs)
        #print formatAllArgs(args, kwargs)
        return func(*args, **kwargs)
    return wrap

def slider(seq, window, start=0):#, stop=None):
    '''assert list(slider([0, 1, 2], 2)) == [ [0,1], [1,2] ]
    assert list(slider('ABCDE', 4)) == [ 'ABCD', 'BCDE' ]
    assert list(slider('ABCDE', 1)) == list('ABCDE')'''
    N = len(seq)
    for idx  in xrange(N-window+1):
        yield seq[idx:idx+window]

filterfst = compose(next, ifilter)
composition = compose_all(sorted, list, slider)

def fromstr(_in):
   lines = filter(str.strip, _in.split('\n'))
   k = int(lines[0])
   s = ''.join(lines[1:])
   return s, k
cfromstr = starcompose(composition, fromstr)
#assert ["AATCC", "ATCCA", "CAATC", "CCAAC", "TCCAA"] == cfromstr(r_in)

#neighobrs = filter(X[:k] == sfx, prefixg)
#NOTE: using generators over lists makes a huge difference.
def make_ovrlp_graph(kmers):
    N = len(kmers)
    ov = len(kmers[0]) - 1
    M = np.zeros((N, N))
    D = dictzip(kmers, xrange(N))
    def update(M, pre, suff):
         sfx, suffg = suff
         nbrs = imap(X[0], ifilter(X[1]==sfx, pre))
         nbr_idxs = map(D.__getitem__, nbrs)
         s_idxs = map(D.__getitem__, suffg)
         M[ s_idxs, nbr_idxs] = 1
    pre, suff = groupby(kmers, X[:ov]), groupby(kmers, X[-ov:])
    pre = izip(kmers, imap(X[:ov], kmers))
    Mupdate = F(update, M, pre)
    #NOTE: starmap drains itertools.groupby
    map(Mupdate, suff)
    return D, M
form = '{0} -> {1}'.format

def printgraph(D, M):
    with open('overlap.txt', 'w') as out:
        idx_kmer_map = dict(map(reversed, D.items()))
        for kmer, idx in D.items():
            for nbr in M[idx].nonzero()[0]:
                print(form(kmer, idx_kmer_map[nbr]), file=out)

getkmers = compose(F(filter, str.strip), mc('split', '\n'))
slv_overlap=compose(starcompose(printgraph, make_ovrlp_graph), getkmers)


exp='''
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG
'''
_in=       '''ATGCG
GCATG
CATGC
AGGCA
GGCAT
'''

def allbut(end, excluding):
    ''' get range(0, end) excluding one element '''
    if excluding < end: return range(start or 0, excluding) + range(excluding + 1, end)
    else: return range(0, end)

def emptygraph(N):
    G = nx.Graph()
    ''' add blank nodes to ease the construction of the graph.'''
    #G.add_nodes_from(zip(xrange(N), [dict(weight=0) for i in xrange(N)]))
    G.add_nodes_from(xrange(N))
    return G


#def prepmatrix(m):
#     M = doublematrix(m)
#     M[np.diag_indices(len(M))] = np.inf
#     return M

def min_matrix_ij(M):
    ''' :return indices (i, j) for the minimum element in the matrix.
    assumes diagonals != 0 or are masked. '''
    return np.unravel_index(M.argmin(), M.shape)

def doublematrix(m):
    n = m.shape[0]
    #M = np.full((n*2, n*2), 0)
    M = np.zeros((n*2, n*2), 0)
    M[0:n, 0:n] = m
    return M

import traceback
def pvar(__x):
        print (traceback.extract_stack(limit=2)[0][3][6:][:-1],"=",__x)













