import networkx as nx
import itertools
from func import compose_all, compose, starcompose, dictzip, pmap
from fn import F, _ as X
import random
from itertools import groupby, starmap, imap, ifilter, izip, ifilterfalse
from operator import methodcaller as mc
import numpy as np
import matplotlib.pyplot as plt
from utils import slider

#G = nx.read_pajek('Sampson.paj')
def drawgraph(G, edgekey='weight', big=False, **kwargs):
    if big: fig = plt.figure(figsize = (15, 10))
    pos=nx.spring_layout(G)
    nx.draw_networkx(G, pos=pos, **kwargs)
    if edgekey:
        edge_labels=dict([((u,v,),d.get(edgekey, ''))
                         for u,v,d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)#, **kwargs)
    plt.show()
#NOTE: requires graphviz
    #nx.write_dot(G,'graph.dot')
#dot -Tpng graph.dot > graph.png

def info_fromkmer(kmer):
        node1, node2 = kmer[:-1], kmer[1:]
        return node1, node2, {'kmer' : kmer}


yield_pathgraph = compose(F(imap, info_fromkmer), slider)
#pathlist = compose(list, yield_pathgraph)
#use reduce
def make_debruijn(s, k=None):
    G = nx.MultiDiGraph()
    if not k:
       G.add_edges_from(imap(info_fromkmer, s))
    else:
        #build straight from string
       G.add_edges_from(yield_pathgraph(s, k))
    return G



'''
set v to some random node.
randomly walk graph until reaching v.
until graphexplored:
    select a node with unexplored edges within the walked cycle.
    create new cycle starting at that edge.
    set the new cycle to be the walked cycle.
'''
filterfalse = compose(list, ifilterfalse)
def walk(G, vstd, cycle, start, current=None, call=0):
    #TODO: I think this leaves out the final step of the cycle.
    if start == current: return vstd, cycle# + tuple([current])
    #NOTE: checking for boolean of 0 is bad here haha
    #_current = start if current else current
    _current = start if current is None else current
    candidates = set(G.edges(_current)) - vstd
    #candidates = filterfalse(vstd.__contains__, G.neighbors(current))
    edge = random.choice(tuple(candidates))
    nn = edge[1]
    return walk(G,  vstd | set([edge]), cycle + tuple([nn]), start, nn, call+1)

filterfst = compose(next, ifilter)
def edges_of_path(G, p):
    return map(X[0]['kmer'], starmap(F(G.get_edge_data), slider(p, 2)))
reconstruct_str = compose_all(''.join, pmap(''.join), edges_of_path)

def e_cycle(G, vstd=set(), cycle=(), call=0):
    ''' find a Eulerian path in a graph by iteratively expanding a cycle.
    requires a mostly-balanced and connected graph.'''
    if len(vstd) == len(G.edges()): return cycle

    def valid(N):
        edges=G.edges(N)
        return not (set(edges) <= vstd)
        #return bool(map(F(filterfalse, vstd.__contains__), edges))
    if not cycle:
        valid_start = random.choice(G.nodes()) # 6
        cycle = tuple([valid_start])
    else:
        valid_start = filterfst(valid, cycle)
    _vstd, _cycle = walk(G,  vstd, cycle, valid_start)
    return e_cycle(G, _vstd, _cycle, call+1)

def k_circular(k):
    kmers = tuple(itertools.product(*(['01']*k)))
    global kg
    kg= make_debruijn(kmers)
    return e_cycle(kg)

k_circular(3)
_in = "TAATGCCATGGGATGTT"
res = make_debruijn(_in, 3)
drawgraph(res, edgekey='kmer')
G = nx.DiGraph()
G.add_edges_from([(0 , 3),
(1 , 0),
(2 , 1),
(2, 6),
(3 , 2),
(4 , 2),
(5 , 4),
(6 , 5),
(6, 8),
(7 , 9),
(8 , 7),
(9 , 6)])

print e_cycle(G)

import subprocess as sp
nx.write_dot(G,'graph.dot')
#sp.check_output("dot -Tpng graph.dot > graph2.png".split())
#sp.check_output("eog graph2.png".split())
