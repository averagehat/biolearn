from trees import get_products, parse_line
from numbers import Number
from func import typecheck, compose, ilen, pmap, dictzip
import networkx as nx
from itertools import starmap, imap, ifilter, takewhile, ifilterfalse
import numpy as np
from testleaves import parsematrix
from fn import _, F
from fn.iters import take, accumulate
from utils import slider
from assembly import drawgraph
from numpy import nan

def to_adj_list(G, edgekey='weight', bothways=True, as_float=False):
    edges = sorted(ifilterfalse(lambda x: x[0] ==x[1], G.edges(data=True)))
    if as_float:
        res = map(lambda x: (x[0], x[1], float(x[-1]['weight'])), edges)
        form = "{0}->{1}:{2:.3f}\n{1}->{0}:{2:.3f}".format
    else:
        res = map(lambda x: (x[0], x[1], int(x[-1]['weight'])), edges)
        form="{0}->{1}:{2}".format if not bothways else "{0}->{1}:{2}\n{1}->{0}:{2}".format
    return starmap(form, res)
adj_str = compose('\n'.join, to_adj_list)

def fst_or_none(func, seq):
    res = filter(func, seq)
    return None if not res else res[0]
filterfst = compose(next, ifilter)

def nondiag(D, i):return list(set(range(len(D))) - set([i]) )
ndiag_perms = compose(get_products, nondiag)

def limbmatch(D, n):
    ''' find nodes i and k such that they satisfie the linear equation:
        D_ik = D_in + D_nk'''
    #print np.isnan(D[n]).all()
    def match(tup):
        i, k = tup
        return D[i,k] == (D[i,n] + D[n,k])
    leaf_perms = ndiag_perms(D, n)
    return filterfst(match, leaf_perms)

@typecheck(np.ndarray, Number)
def limbmin(D, j):
    '''Same as limblen but also return i & k; not useful.'''
    def LL(i, k): return (D[i,j] + D[j,k] - D[i,k])/2, i, k
    nondiag = set(range(len(D))) - set([j])
    leaf_perms = get_products(nondiag)
    dst, i, k =  min(starmap(LL, leaf_perms), key=_[0])
    return dst, i, k

@typecheck(np.ndarray, Number)
def limblen(D, j):
    '''given a distance matrix *D* and leaf index *j*, return the
    length of the limb which *j* could be added to the simple tree
    represented by *D*.'''
    #NOTE: I had a problem with this function returning NAN.
    def LL(i, k): return (D[i,j] + D[j,k] - D[i,k])/2
    nondiag = set(range(len(D))) - set([j])
    leaf_perms = get_products(nondiag)
    return min(ifilterfalse(np.isnan, starmap(LL, leaf_perms)))



def simplepath(G,  dest, vstd=set(), path=(), source=0):
    ''' find a (unique) path in a simple graph. '''
    if source == dest:
        return path + tuple([dest])
    else:
        nextpath = F(simplepath, G, dest, vstd | set([source]), path + tuple([source]))
        nbrs = set(G.neighbors(source)) - vstd
        results = filter(bool, map(nextpath, nbrs))
        return None if not results else results[0]


def nd_matching_path_wt(G, path, dist):
    '''given a path, and a graph, find the node which is *dist* distance along the path.
    if no such path exists, return the surrounding nodes. '''
    pathinfo = starmap(G.get_edge_data, slider(path, 2))
    weights = map(_['weight'], pathinfo)
    acc_weights =  list(accumulate(weights))
    if dist not in acc_weights:
        #then it is between two nodes
        idx=ilen(takewhile(_<dist, acc_weights))
        return path[idx], path[idx+1]
    return path[acc_weights.index(dist)]

def node_along_path_matching_weight(G, start, end, dist):
    ''' see nd_matching_path_wt'''
    path = simplepath(G, end, source=start)
    return nd_matching_path_wt(G, path, dist)

def np_data(M):
    ''' return a networkx friendly list of tuples (key, val, weight)'''
    a = []
    for i in xrange(len(M)):
        for j in xrange(M.shape[1]):
            if not np.isnan(M[i,j]):
                a.append((i, j, int(M[i, j])))
    return a
#rans = (i-1 for i in [1+1, 3+1, 2+1, 0])
get_weight = compose(_['weight'], nx.DiGraph.get_edge_data)
#NOTE: nx constructor does not handle np.masked array properly
def add_phylo(D, n, draw=False):
    '''build a simple unrooted phylogenetic tree from a distance matrix.'''
    if n == 1: return fromnp(D) #nx.from_numpy_matrix(D)
    '''compute the length of the limb we will attach n to.'''
    ll = limblen(D, n)
    D[n] -= ll; D[:, n] -= ll
    _i, _k = limbmatch(D, n)
    x = D[_i, n]
    '''remove n from the distance matrix.'''
    D[n] = D[:, n] = nan
    #recursive tree creation; D.copy unecessary, used for debugging
    T = add_phylo(D.copy(), n-1)
    ''' find the nodes to put the new internal node (witch we attach to n) between. '''
    v = node_along_path_matching_weight(T, _i, _k, x)
    if len(v) == 1:
    #NOTE: this never happens
         T.add_edge(v, n, {'weight' : ll})
    else:
         ''' add node N and any needed internal node to the graph.'''
         if draw: drawgraph(T)
         #i and k are the nodes that the inserted node lies between.
         ''' there are two cases here. '''
         #TODO: are they necessary?
         i, k = v
         if _i == i:
             z = x
         else:
             #NOTE: shortest_path_length will default to weight=1 per edge if not given the weight kwarg.
             dii = nx.shortest_path_length(T, _i, i, weight='weight')
             z = x - dii
         klen = get_weight(T, i, k) - z
#         print '_i: %s, _k: %s, i: %s, k: %s' % (_i, _k, i, k)
#         print "ll: %s , x: %s , z: %s , klen: %s , dii: %s" % (ll, x, z, klen, dii)
         nv = D.shape[0] + n - 2
         if T.has_edge(i, k):
             T.remove_edge(i, k)
             T.remove_edge(k, i)
         #TODO: factor out
         ''' add the new internal node and update the relevant lengths.'''
         T.add_edge(n, nv, {'weight' : ll})
         T.add_edge(nv, n, {'weight' : ll})
         T.add_edge(i, nv, {'weight' : z})
         T.add_edge(nv, i, {'weight' : z})
         T.add_edge(k, nv, {'weight' : klen})
         T.add_edge(nv, k, {'weight' : klen})
    return T
#compute limblength #and i, k using limbmin #substract limblength from D[n] #set i, k to be limb_match(n) #set ex to be distance between i and n (D[i, n]) #mask D[n] by setting nan.  #T = recur(D, n) # get v out of T

'''
v may not be a new node. if it is not a new node, it cannot become an internal
node, because it represents a present species.
if v is not a new node, it may need to be replaced and re-connected with a new internal node.
'''

def fixed_lines(raw): return ifilter(str.strip, raw.splitlines())
def make_nxgraph(data, edgekey='weight'):
    ''':param data: list of tuples of form (from, to, edge_info)
        :return nx.DiGraph with *data* edges.'''
    G = nx.DiGraph()
    _data = [(d[0], d[1], {edgekey : d[2]}) for d in data]
    G.add_edges_from(_data)
    return G
parse_lines = compose(F(imap, parse_line), fixed_lines)
G_from_adj_list = compose(make_nxgraph, parse_lines)

fromnp = compose(make_nxgraph, np_data)
'''
limb length, given a distance matrix, finds the hypothetical length to what would
be the closest node, wether or not this node exists in the distance matrix.
It does this by applying the equation LL to all other possible leaves.
'''

def test_basic():
    dt=[('weight',float),('cost',int)]
    _in1 = np.array([[0, 21, 22],
                   [21, 0, 13],
                   [22, 13, 0]], dtype=float)

    _in = np.matrix([[0 ,  13,  21,  22],
[13,  0 ,  12,  13],
[21,  12,  0 ,  13],
[22,  13,  13,  0]], dtype=float)#dtype=[('weight', float), ('cost', int)])
    res2 = add_phylo(_in, 3)
    drawgraph(res2)
    #drawgraph(res2)
    #drawgraph(add_phylo(_in1, 2))
    M = np.matrix([[0, 11, 10, 9, 15],
     [11, 0, 3, 12, 18],
     [10, 3, 0, 11, 17],
     [9, 12, 11, 0, 8],
     [15, 18, 17, 8, 0]], dtype=float)
    global big
    big = add_phylo(M, len(M) - 1)
    drawgraph(big)




M= np.array([[   0 ,  13,  21,  22],
[   13,  0 ,  12,  13],
[   21,  12,  0 ,  13],
[   22,  13,  13,  0]])
j = 1
assert limblen(M, j) == 2
ac = test_basic()
ex = add_phylo(np.array([[0, 21], [21, 0]]), 1)
pathed = simplepath(ex, 1)
def parseM(raw):
    '''parse & return a space-seperated matrix.'''
    _in = filter(bool, raw.split('\n'))
    return  np.matrix(map(pmap(float), map(str.split, _in)))


M= parseM('''0 3036 4777 1541 2766 6656 2401 4119 7488 4929 5344 3516 1485 6392 2066 3216 7008 7206 1187 6491 3379 6262 6153 4927 6670 4997 9010 5793 9032
3036 0 6323 3087 4312 8202 1619 5665 9034 2205 6890 966 3031 7938 3612 492 4284 8752 2023 8037 4925 3538 3429 6473 3946 2273 10556 7339 10578
4777 6323 0 4054 2571 3455 5688 1972 4287 8216 2143 6803 4126 3191 3183 6503 10295 4005 4474 3290 2964 9549 9440 1726 9957 8284 5809 2592 5831
1541 3087 4054 0 2043 5933 2452 3396 6765 4980 4621 3567 890 5669 1343 3267 7059 6483 1238 5768 2656 6313 6204 4204 6721 5048 8287 5070 8309
2766 4312 2571 2043 0 4450 3677 1913 5282 6205 3138 4792 2115 4186 1172 4492 8284 5000 2463 4285 1173 7538 7429 2721 7946 6273 6804 3587 6826
6656 8202 3455 5933 4450 0 7567 3851 1082 10095 1872 8682 6005 1992 5062 8382 12174 800 6353 1047 4843 11428 11319 3053 11836 10163 2604 2545 2626
2401 1619 5688 2452 3677 7567 0 5030 8399 3512 6255 2099 2396 7303 2977 1799 5591 8117 1388 7402 4290 4845 4736 5838 5253 3580 9921 6704 9943
4119 5665 1972 3396 1913 3851 5030 0 4683 7558 2539 6145 3468 3587 2525 5845 9637 4401 3816 3686 2306 8891 8782 2122 9299 7626 6205 2988 6227
7488 9034 4287 6765 5282 1082 8399 4683 0 10927 2704 9514 6837 2824 5894 9214 13006 1172 7185 1879 5675 12260 12151 3885 12668 10995 2144 3377 2166
4929 2205 8216 4980 6205 10095 3512 7558 10927 0 8783 2859 4924 9831 5505 1891 3719 10645 3916 9930 6818 2973 2864 8366 3381 1708 12449 9232 12471
5344 6890 2143 4621 3138 1872 6255 2539 2704 8783 0 7370 4693 1608 3750 7070 10862 2422 5041 1707 3531 10116 10007 1741 10524 8851 4226 1233 4248
3516 966 6803 3567 4792 8682 2099 6145 9514 2859 7370 0 3511 8418 4092 1146 4938 9232 2503 8517 5405 4192 4083 6953 4600 2927 11036 7819 11058
1485 3031 4126 890 2115 6005 2396 3468 6837 4924 4693 3511 0 5741 1415 3211 7003 6555 1182 5840 2728 6257 6148 4276 6665 4992 8359 5142 8381
6392 7938 3191 5669 4186 1992 7303 3587 2824 9831 1608 8418 5741 0 4798 8118 11910 2542 6089 1827 4579 11164 11055 2789 11572 9899 4346 2281 4368
2066 3612 3183 1343 1172 5062 2977 2525 5894 5505 3750 4092 1415 4798 0 3792 7584 5612 1763 4897 1785 6838 6729 3333 7246 5573 7416 4199 7438
3216 492 6503 3267 4492 8382 1799 5845 9214 1891 7070 1146 3211 8118 3792 0 3970 8932 2203 8217 5105 3224 3115 6653 3632 1959 10736 7519 10758
7008 4284 10295 7059 8284 12174 5591 9637 13006 3719 10862 4938 7003 11910 7584 3970 0 12724 5995 12009 8897 1442 2699 10445 1088 3447 14528 11311 14550
7206 8752 4005 6483 5000 800 8117 4401 1172 10645 2422 9232 6555 2542 5612 8932 12724 0 6903 1597 5393 11978 11869 3603 12386 10713 2694 3095 2716
1187 2023 4474 1238 2463 6353 1388 3816 7185 3916 5041 2503 1182 6089 1763 2203 5995 6903 0 6188 3076 5249 5140 4624 5657 3984 8707 5490 8729
6491 8037 3290 5768 4285 1047 7402 3686 1879 9930 1707 8517 5840 1827 4897 8217 12009 1597 6188 0 4678 11263 11154 2888 11671 9998 3401 2380 3423
3379 4925 2964 2656 1173 4843 4290 2306 5675 6818 3531 5405 2728 4579 1785 5105 8897 5393 3076 4678 0 8151 8042 3114 8559 6886 7197 3980 7219
6262 3538 9549 6313 7538 11428 4845 8891 12260 2973 10116 4192 6257 11164 6838 3224 1442 11978 5249 11263 8151 0 1953 9699 1104 2701 13782 10565 13804
6153 3429 9440 6204 7429 11319 4736 8782 12151 2864 10007 4083 6148 11055 6729 3115 2699 11869 5140 11154 8042 1953 0 9590 2361 2592 13673 10456 13695
4927 6473 1726 4204 2721 3053 5838 2122 3885 8366 1741 6953 4276 2789 3333 6653 10445 3603 4624 2888 3114 9699 9590 0 10107 8434 5407 2190 5429
6670 3946 9957 6721 7946 11836 5253 9299 12668 3381 10524 4600 6665 11572 7246 3632 1088 12386 5657 11671 8559 1104 2361 10107 0 3109 14190 10973 14212
4997 2273 8284 5048 6273 10163 3580 7626 10995 1708 8851 2927 4992 9899 5573 1959 3447 10713 3984 9998 6886 2701 2592 8434 3109 0 12517 9300 12539
9010 10556 5809 8287 6804 2604 9921 6205 2144 12449 4226 11036 8359 4346 7416 10736 14528 2694 8707 3401 7197 13782 13673 5407 14190 12517 0 4899 1758
5793 7339 2592 5070 3587 2545 6704 2988 3377 9232 1233 7819 5142 2281 4199 7519 11311 3095 5490 2380 3980 10565 10456 2190 10973 9300 4899 0 4921
9032 10578 5831 8309 6826 2626 9943 6227 2166 12471 4248 11058 8381 4368 7438 10758 14550 2716 8729 3423 7219 13804 13695 5429 14212 12539 1758 4921 0'''  )

G1 = add_phylo(M, len(M) - 1)
print adj_str(G1)
M=parseM('''
0 5556 2330 5253 2651 7603 3471 5720 5121 966 6079 7178 2413 1925 2469 2962 4514 4490 1737 3505 4380 6552 7002 3064 6629 2432
5556 0 7512 1541 7833 12785 2681 10902 1995 6148 11261 12360 3787 5303 7651 3944 9696 1960 6919 8687 9562 11734 12184 8246 11811 3974
2330 7512 0 7209 2205 7157 5427 5274 7077 2452 5633 6732 4369 3881 2023 4918 4068 6446 2013 3059 3934 6106 6556 2618 6183 4388
5253 1541 7209 0 7530 12482 2378 10599 1692 5845 10958 12057 3484 5000 7348 3641 9393 1657 6616 8384 9259 11431 11881 7943 11508 3671
2651 7833 2205 7530 0 6922 5748 5039 7398 2773 5398 6497 4690 4202 1788 5239 3833 6767 2334 2824 3699 5871 6321 2383 5948 4709
7603 12785 7157 12482 6922 0 10700 3129 12350 7725 3266 1033 9642 9154 6066 10191 3233 11719 7286 4426 3899 1393 1989 5627 1904 9661
3471 2681 5427 2378 5748 10700 0 8817 2246 4063 9176 10275 1702 3218 5566 1859 7611 1615 4834 6602 7477 9649 10099 6161 9726 1889
5720 10902 5274 10599 5039 3129 8817 0 10467 5842 1605 2704 7759 7271 4183 8308 1350 9836 5403 2543 2016 2078 2528 3744 2155 7778
5121 1995 7077 1692 7398 12350 2246 10467 0 5713 10826 11925 3352 4868 7216 3509 9261 1525 6484 8252 9127 11299 11749 7811 11376 3539
966 6148 2452 5845 2773 7725 4063 5842 5713 0 6201 7300 3005 2517 2591 3554 4636 5082 1859 3627 4502 6674 7124 3186 6751 3024
6079 11261 5633 10958 5398 3266 9176 1605 10826 6201 0 2841 8118 7630 4542 8667 1709 10195 5762 2902 2375 2215 2665 4103 2292 8137
7178 12360 6732 12057 6497 1033 10275 2704 11925 7300 2841 0 9217 8729 5641 9766 2808 11294 6861 4001 3474 968 1564 5202 1479 9236
2413 3787 4369 3484 4690 9642 1702 7759 3352 3005 8118 9217 0 2160 4508 1193 6553 2721 3776 5544 6419 8591 9041 5103 8668 831
1925 5303 3881 5000 4202 9154 3218 7271 4868 2517 7630 8729 2160 0 4020 2709 6065 4237 3288 5056 5931 8103 8553 4615 8180 2179
2469 7651 2023 7348 1788 6066 5566 4183 7216 2591 4542 5641 4508 4020 0 5057 2977 6585 2152 1968 2843 5015 5465 1527 5092 4527
2962 3944 4918 3641 5239 10191 1859 8308 3509 3554 8667 9766 1193 2709 5057 0 7102 2878 4325 6093 6968 9140 9590 5652 9217 1380
4514 9696 4068 9393 3833 3233 7611 1350 9261 4636 1709 2808 6553 6065 2977 7102 0 8630 4197 1337 810 2182 2632 2538 2259 6572
4490 1960 6446 1657 6767 11719 1615 9836 1525 5082 10195 11294 2721 4237 6585 2878 8630 0 5853 7621 8496 10668 11118 7180 10745 2908
1737 6919 2013 6616 2334 7286 4834 5403 6484 1859 5762 6861 3776 3288 2152 4325 4197 5853 0 3188 4063 6235 6685 2747 6312 3795
3505 8687 3059 8384 2824 4426 6602 2543 8252 3627 2902 4001 5544 5056 1968 6093 1337 7621 3188 0 1203 3375 3825 1529 3452 5563
4380 9562 3934 9259 3699 3899 7477 2016 9127 4502 2375 3474 6419 5931 2843 6968 810 8496 4063 1203 0 2848 3298 2404 2925 6438
6552 11734 6106 11431 5871 1393 9649 2078 11299 6674 2215 968 8591 8103 5015 9140 2182 10668 6235 3375 2848 0 938 4576 853 8610
7002 12184 6556 11881 6321 1989 10099 2528 11749 7124 2665 1564 9041 8553 5465 9590 2632 11118 6685 3825 3298 938 0 5026 1303 9060
3064 8246 2618 7943 2383 5627 6161 3744 7811 3186 4103 5202 5103 4615 1527 5652 2538 7180 2747 1529 2404 4576 5026 0 4653 5122
6629 11811 6183 11508 5948 1904 9726 2155 11376 6751 2292 1479 8668 8180 5092 9217 2259 10745 6312 3452 2925 853 1303 4653 0 8687
2432 3974 4388 3671 4709 9661 1889 7778 3539 3024 8137 9236 831 2179 4527 1380 6572 2908 3795 5563 6438 8610 9060 5122 8687 0'''   )

G2 = add_phylo(M, len(M) - 1)
print adj_str(G2)
exp = open('expaddphylo.txt').read().strip()
assert adj_str(G2, bothways=False) == exp
'Test Passed.'

names='Cow Pig Horse   Mouse   Dog Cat Turkey  Civet   Human'.split()
mapping = dict(enumerate(names))


real = parseM(''' 0   295 306 497 1081    1091    1003    956 954
 295 0   309 500 1084    1094    1006    959 957
   306 309     0   489 1073    1083    995 948 946
   497 500 489 0   1092    1102    1014    967 965
 1081    1084    1073    1092    0   818 1056    1053    1051
 1091    1094    1083    1102    818 0   1066    1063    1061
1003    1006    995 1014    1056    1066    0   975 973
956 959 948 967 1053    1063    975 0   16
954 957 946 965 1051    1061    973 16  0''')


RG = add_phylo(real, len(real) - 1)
RG = nx.relabel_nodes(RG, mapping)
