from __future__ import print_function
from func import compose_all, compose, starcompose, dictzip, ilen
from fn import F, _ as X
from itertools import groupby, starmap, imap, ifilter, izip
from operator import methodcaller as mc, ne
import numpy as np

def quantify(iterable, pred=bool):
    '''https://docs.python.org/2/library/itertools.html#recipes
    "Count how many times the predicate is true"'''
    return sum(imap(pred, iterable))

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
    N = len(seq)
    for idx  in xrange(N-window+1):
        yield seq[idx:idx+window]

filterfst = compose(next, ifilter)

assert list(slider([0, 1, 2], 2)) == [ [0,1], [1,2] ]

assert list(slider('ABCDE', 4)) == [ 'ABCD', 'BCDE' ]


assert list(slider('ABCDE', 1)) == list('ABCDE')
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

#slv_overlap(open('../../Downloads/rosalind_4b.txt').read())

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

