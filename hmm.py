''' hidden markov models'''
import numpy as np
import string
import pandas as pd
from func import starcompose, partial2
import math
import StringIO
from utils import slider
from fn import F, _
import operator
from fn.iters import filterfalse, map
d = {'A' : 0, 'B' : 1}
ed = {'x':0, 'y':1,'z':2}
X = _
def pr_hidden_path(transition_matrix, hiddenpath):
    ''' What is the probability, given a transition matrix between hidden states,
    That a series of hidden states (a hidden path) could arise?'''
    M, hp = transition_matrix, hiddenpath

    def acc_transition(acc, trans):
        pr = M[ d[trans[0]], d[trans[1]] ]
        return pr * acc
    ''' Reduction by product, where the input is the transition (reflected by slider 2) '''
    return reduce(acc_transition, slider(hp, 2), 0.5)



a = np.array([[ 0.377,   0.623], [0.26,    0.74]])
p = 'ABABBBAAAA'


assert (0.000384928691755 -   pr_hidden_path(a, p))  < 0.000001


p = 'BBABBBABBAABABABBBAABBBBAAABABABAAAABBBBBAABBABABB'
a= np.array([[     0.863 ,  0.137   ], [     0.511 ,  0.489]])

assert (3.26233331904e-21 -   pr_hidden_path(a, p)  ) < 0.000001
p = 'BBABAABABABABABBAAABABBBAAAAAABBBBBBBABABBAABAABAB'

a = np.array([[0.339,    0.661   ], [  0.64,    0.36]])

def fixm(df):
    return df.set_index(df[df.columns[0]]).drop(df.columns[0], 1)
#df['x']['A']
readtabs  = F(pd.read_csv, delimiter='\t')
readM  = F(fixm) << readtabs << F(StringIO.StringIO)
filterfalseL = F(list) << filterfalse
splitter = '--------'
groups = F(filter, str.split) << F(filter, X != '--------\n')
getgroups = groups <<  open
def parsexi(fn):
    groups = getgroups(fn)
    numhdnstates = len(groups[4].split(' ')) + 2
    rawm =  ''.join(groups[-numhdnstates:])
    emission_matrix = readM(rawm)
    expressions, hp= groups[0], groups[2]
    return expressions.strip(), hp.strip(), emission_matrix
fn = '/tmp/dataset_11594_2.txt'
#print parsexi(fn)

def _pr_outcome(emissions, hp, EM):
   # emissions_pr = F(map, EM[X][X])
    product = F(reduce, operator.mul)
   # outcome_pr = product << emissions_pr
    outcome_pr = product << F(map, lambda x, y: EM[x][y])

    return outcome_pr(emissions, hp)

prfull = starcompose(_pr_outcome, parsexi)
AE = np.testing.assert_almost_equal
AE(prfull(fn), 3.59748954746e-06)

AE(prfull('datas/foo'), 3.42316482177e-35)


d = {'A' : 0, 'B' : 1}
od = {0 : 'A', 1 : 'B', 2 : 'C', 3 : 'D'  }
ed = {'x':0, 'y':1,'z':2}
TM = np.array([[0.641,     0.359 ],
               [0.729,     0.271 ]])
EM = np.array([[0.117,     0.691,   0.192],
               [0.097,     0.42 ,   0.483]])
ems = 'xyxzzxyxyy'
tbl = string.maketrans('xyz', '012')
toints = F(np.array) << list << F(map, int) << F(X.call('translate', tbl))

def max_argmax(A):
    agmx = A.argmax()
    return A[agmx], agmx

def _viterbi_log(emissions, TM, EM, scores, n=0):
    if n == len(emissions):
        return scores
    ems = EM[:, emissions[n]]
    if n == 0:
        scores[:, n] = np.log2(0.5 * ems)
        return _viterbi(emissions, TM, EM, scores, n+1)
    #import ipdb; ipdb.set_trace()
    candidates = np.log2(TM) + scores[:, n-1]
    scores[:, n] = np.log2(ems[:, None])+ candidates.max(axis=1)
    return _viterbi(emissions, TM, EM, scores, n+1)

def viterbi(emissions, TM, EM):
    ems = toints(emissions)
    scores = np.full( (TM.shape[0], len(ems)), -np.inf )
    return _viterbi(ems, TM, EM, scores, 0)

#    def accviterbi(acc, em):
#        scores, n = acc
#        candidates = scores[n-1] * TM * ems
#        scores[n] = candidates.max(axis=1)
#        return _viterbi(emissions, TM, EM, scores, n+1)
#
#
#    pass
#log2 = partial2(math.log, 2)



def _viterbi(emissions, TM, EM, scores, n=0):
    if n == len(emissions):
        return scores
    ems = EM[:, emissions[n]]
    #import ipdb; ipdb.set_trace()
    if n == 0:
        #scores[:, n] = 0.5 * ems
        #start with equal probability at each node
        scores[:, n] = 1/float(len(TM)) * ems
        return _viterbi(emissions, TM, EM, scores, n+1)

    #candidates = scores[:, n-1] * TM.T * ems
    scores[:, n] = (scores[:, n-1, None] * TM * ems).max(axis=0)
    #scores[:, n] = candidates.max(axis=1)
    return _viterbi(emissions, TM, EM, scores, n+1)















def backtrk(scores, TM):
    picked = scores[:, -1].argmax()
    res = [picked]
    for i in xrange(2, scores.shape[1]+1):
        picked = (scores[:, -i] * TM[:, picked]).argmax()
        res.append(picked)
    return list(reversed(res))




score =  viterbi(ems, TM, EM)
res = backtrk(score, TM)
final= ''.join(map(od.__getitem__, res))
assert final == 'AAABBAAAAA'



s = 'zxxxxyzzxyxyxyzxzzxzzzyzzxxxzxxyyyzxyxzyxyxyzyyyyzzyyyyzzxzxzyzzzzyxzxxxyxxxxyyzyyzyyyxzzzzyzxyzzyyy'
tm= np.array([[  0.634,   0.366],   [  0.387,   0.613],   ])
em = np.array([[  0.532,   0.226,   0.241],   [  0.457,   0.192,   0.351]])


score =  viterbi(s , tm, em)
res = backtrk(score, tm)
final =   ''.join(map(od.__getitem__, res))
assert final == 'AAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBAAA'



s = 'xxyxxzzzzxyxyxzxyzyxzyzzzzzyyyxyyxxxzyxxxzzzxzxzyzxxzyyyzxxzzxzyzzxzzyzzzzxyyyzyxyzzyxzzxxyxxzyyzyzx'

tm = np.array([[0.018,   0.381,   0.16  ,  0.441   ],
[0.361,   0.346,   0.162 ,  0.131   ],
[0.319,   0.135,   0.451 ,  0.095   ],
[0.108,   0.442,   0.075 ,  0.375   ]])

em = np.array([ [0.298,   0.324,   0.378 ,          ],
[0.261,   0.375,   0.364 ,          ],
[0.195,   0.593,   0.212 ,          ],
[0.286,   0.184,   0.53             ]])

s = 'yxxxyyxzzxzyxyyyyzxyyyzzzxzyxyxzzyyzxzzyzyxyyyzzyzyxzzzyzxyyzzyxxyyxyyyxzyyxyyxyyyzzxxzxyxxxyxyyyzxy'
tm=np.array   ([[0.348 ,  0.215 ,  0.276,   0.161   ],
[0.742 ,  0.027 ,  0.188 ,  0.043   ],
[0.347 ,  0.348 ,  0.146 ,  0.159   ],
[0.016 ,  0.366 ,  0.258 ,  0.36    ]])

em = np.array([[0.542 ,  0.256 ,  0.202   ],
[0.474 ,  0.5   ,   0.026   ],
[0.257 ,  0.294 ,  0.449   ],
[0.265 ,  0.117 ,  0.618]])
score =  viterbi(s , tm, em)
res = backtrk(score, tm)
final =   ''.join(map(od.__getitem__, res))





_tbl = string.maketrans('ACGT', '0123')
_toints = F(np.array) << list << F(map, int) << F(X.call('translate', _tbl))

def forward(emissions, TM, EM, tfunc=toints):
    ems = tfunc(emissions)
    scores = np.zeros( (TM.shape[0], len(ems)))
    return _forward(ems, TM, EM, scores, 0)

def _forward(emissions, TM, EM, scores, n=0):
    if n == len(emissions):
        return scores
    ems = EM[:, emissions[n]]
    #import ipdb; ipdb.set_trace()
    if n == 0:
        #start with equal probability at each node
        scores[:, n] = 1/float(len(TM)) * ems
        return _forward(emissions, TM, EM, scores, n+1)
    scores[:, n] = (scores[:, n-1, None] * TM * ems).sum(axis=0)
#    for i in xrange(scores.shape[-1]):
#        #print TM[i], em[i], scores[i, n]
#        scores[i, n] = (TM[i] * em[i] + scores[i, n]).sum()
    return _forward(emissions, TM, EM, scores, n+1)

A = np.array
s = 'GGCA'
tm = A( [ [0.5, .5],
          [.4, .6]])
em = A([[ 0.2, .3, .3, .2],
       [ 0.3, .2, .2, .3]])
score =  forward(s , tm, em, _toints)
#res = backtrk(score, tm)
#print score
s = 'xzyyzzyzyy'
tm= A([[0.303,   0.697], [0.831,   0.169], ])
em = A([[0.533,   0.065 ,  0.402 ],
[0.342,   0.334 ,  0.324 ]])
score =  forward(s , tm, em)
#print score[:, -1].sum()
#res = backtrk(score, tm)

s = 'yyzzzzxzxxxxzxzzzxzzyyzyzyzxyzzyyyzzzzzyxyyyzzzyzxzzzxyxyyxzyzxzzyzxxzyxxzxzyxzxxyzyyxzzzyzzxyyyzzyz'
tm = A([[0.6  ,  0.055 ,  0.345   ],
[0.248 ,  0.275 ,  0.477   ],
[0.327 ,  0.14  ,  0.533   ],])
em = A([[0.199 ,  0.308 ,  0.493   ],
[0.19  ,  0.392 ,  0.418   ],
[0.512 ,  0.155 ,  0.333],])
score =  forward(s , tm, em)
print score[:, -1].sum()
s = 'zxxxzyyxyzyxyyxzzxzyyxzzxyxxzyzzyzyzzyxxyzxxzyxxzxxyzzzzzzzxyzyxzzyxzzyzxyyyyyxzzzyzxxyyyzxyyxyzyyxz'
tm = A([[0.994,   0.006]   , [0.563,   0.437  ]])
em = A([[0.55,    0.276 ,  0.174   ], [0.311,   0.368,   0.321]])
score =  forward(s , tm, em)
res =  score[:, -1].sum()
AE(res, 4.08210708381e-55)

s = 'xzxzxzxzxzzyxyyyyxzyyxzyxyzxzyxzxxzyyxyxxyxzyxzyzxxzxzyyxzxxzzyzzzxyzxzzyyyxyzyzxzxyzzyxxyxyzyzyyzyz'
tm = A([[0.503 ,  0.18  ,  0.317,],
[0.326 ,  0.314 ,  0.36 ,],
[0.331 ,  0.163 ,  0.506,]]   )
em = A([ [0.308 ,  0.09  ,  0.602,],
[0.658 ,  0.164 ,  0.178,],
[0.08  ,  0.552 ,  0.368,]])
score =  forward(s , tm, em)
res =  score[:, -1].sum()




s = 'xzxzyyzyzxzyxyxyxxxxxxyzxxxzzzzzxzzxxyxzyxxyxyxzyxzxzyxyzxxzyyxzzyzzxzzzyzzyzyxxxzzyzzxzyyxyyxyxxxxz'
tm = A([[0.34 ,   0.66   ], [0.374,   0.626  ], ])
em = A([[0.454,   0.235  , 0.311   ], [0.187,   0.279  , 0.534]])
score =  forward(s , tm, em)
res =  score[:, -1].sum()
