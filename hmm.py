''' hidden markov models'''
from __future__ import division
from fplearn import Ratio
import numpy as np
import string
from itertools import dropwhile, starmap
import pandas as pd
from func import starcompose, partial2
import math
import StringIO
from utils import slider
from fn import F, _
import operator
from fn.iters import filterfalse

d = {'A' : 0, 'B' : 1}
od = {0 : 'A', 1 : 'B', 2 : 'C', 3 : 'D'  }
ed = {'x':0, 'y':1,'z':2}
X = _
def fixm(df):
    return df.set_index(df[df.columns[0]]).drop(df.columns[0], 1)

readtabs  = F(pd.read_csv, delimiter='\t')
readM  = F(fixm) << readtabs << F(StringIO.StringIO)
filterfalseL = F(list) << filterfalse
splitter = '--------'
groups = F(filter, str.split) << F(filter, X != '--------\n')
getgroups = groups <<  open

def pr_hidden_path(transition_matrix, hiddenpath):
    ''' What is the probability, given a transition matrix between hidden states,
    That a series of hidden states (a hidden path) could arise?'''
    M, hp = transition_matrix, hiddenpath

    def acc_transition(acc, trans):
        pr = M[ d[trans[0]], d[trans[1]] ]
        return pr * acc
    ''' Reduction by product, where the input is the transition (reflected by slider 2) '''
    return reduce(acc_transition, slider(hp, 2), 0.5)

def parsexi(fn):
    groups = getgroups(fn)
    numhdnstates = len(groups[4].split(' ')) + 2
    rawm =  ''.join(groups[-numhdnstates:])
    emission_matrix = readM(rawm)
    expressions, hp= groups[0], groups[2]
    return expressions.strip(), hp.strip(), emission_matrix

def _pr_outcome(emissions, hp, EM):
    product = F(reduce, operator.mul)
    outcome_pr = product << F(map, lambda x, y: EM[x][y])
    return outcome_pr(emissions, hp)

prfull = starcompose(_pr_outcome, parsexi)
tbl = string.maketrans('xyz', '012')
toints = F(np.array) << list << F(map, int) << F(X.call('translate', tbl))

def viterbi(emissions, TM, EM):
    ems = toints(emissions)
    scores = np.full( (TM.shape[0], len(ems)), -np.inf )
    return _viterbi(ems, TM, EM, scores, 0)

#    def accviterbi(acc, em):
#        scores, n = acc
#        candidates = scores[n-1] * TM * ems
#        scores[n] = candidates.max(axis=1)
#        return _viterbi(emissions, TM, EM, scores, n+1)



def _viterbi(emissions, TM, EM, scores, n=0):
    if n == len(emissions):
        return scores
    ems = EM[:, emissions[n]]
    if n == 0:
        #start with equal probability at each node
        scores[:, n] = 1/float(len(TM)) * ems
        return _viterbi(emissions, TM, EM, scores, n+1)

    #candidates = scores[:, n-1] * TM.T * ems
    scores[:, n] = (scores[:, n-1, None] * TM * ems).max(axis=0)
    return _viterbi(emissions, TM, EM, scores, n+1)

def backtrk(scores, TM):
    picked = scores[:, -1].argmax()
    res = [picked]
    for i in xrange(2, scores.shape[1]+1):
        picked = (scores[:, -i] * TM[:, picked]).argmax()
        res.append(picked)
    return list(reversed(res))

_tbl = string.maketrans('ACGT', '0123')
_toints = F(np.array) << list << F(map, int) << F(X.call('translate', _tbl))

def forward(emissions, TM, EM, tfunc=toints):
    '''Find the probability that a series of tokens (emissions) matches
    a given model (described by transition and emission matrices)'''
    ''' use the sum where viterbi uses max!'''
    ems = tfunc(emissions)
    scores = np.zeros( (TM.shape[0], len(ems)))
    return _forward(ems, TM, EM, scores, 0)

def _forward(emissions, TM, EM, scores, n=0):
    if n == len(emissions):
        return scores
    ems = EM[:, emissions[n]]
    if n == 0:
        #start with equal probability at each node
        scores[:, n] = 1/float(len(TM)) * ems
        return _forward(emissions, TM, EM, scores, n+1)
    ''' use the sum where viterbi uses max!'''
    scores[:, n] = (scores[:, n-1, None] * TM * ems).sum(axis=0)
    return _forward(emissions, TM, EM, scores, n+1)


from fn.iters import range
import itertools
EMPTY = '-'
HSTATES = ('I', 'M', 'D')
mapjoin = F(map, F(''.join) << F(map, str))
make_hstates = mapjoin << F(itertools.product, HSTATES)  << range
#makedf = F(lambda x, i: pd.DataFrame(i, index=x, columns=x)) << F(list) << make_hstates
_makedf = F(lambda x: pd.DataFrame(0, index=['S'] + x + ['E'], columns= ['S'] + x  + ['E']))
st = lambda x, y : '%s%s' % (x, y)
def fixperms(perms):
    print len(perms)
    n = len(perms)//4
    bads = ['M0', 'D0'] #+ map(partial2(st, n), ['M', 'I', 'D'])
    print bads
    return list( set(perms) - set(bads) )

makedf = _makedf << fixperms << F(list) << make_hstates
#makedf =_makedf << F(list) << make_hstates

def make_hmm(C, theta):
    N = C.shape[0]
    empties = C == EMPTY
    nempties = C != EMPTY
    is_insert  = ( (C == EMPTY).sum(axis=0) / N ) > theta
    #deletion = empties.T * ~is_insert
    #match = nempties.T * ~is_insert
    #insert = nempties.T * is_insert
    deletion = (empties * ~is_insert).T
    match =    (nempties * ~is_insert).T
    insert =   (nempties * is_insert).T
    n_deletion, n_match, n_insert = map(F(np.sum, axis=1), [deletion, match, insert])
    d = { "M" : n_match, "D" : n_deletion, "I" : n_insert}
    d2 = { "M" : match, "D" : deletion, "I" : insert}
    cells = C.shape[1]
    df = makedf(cells, 0)
    #import ipdb; ipdb.set_trace()
    print df
    c=0
    for i in range(cells-1):
        #TODO: this doesn't work because M0 -> M2 not M0 -> M1; and this looks forward only one ahead.
        for name1, M in d.items():
            selected = set()
            c+=1
            for name2, M2 in d.items():
                #df['%s%s' %(name1, i), '%s%s' % (name2, i+1)] = M[i] + M2[i+1]
                #if name1 == 'M' and i == 0: import ipdb; ipdb.set_trace()
                i1, i2 = '{}{}'.format(name1, i), '{}{}'.format(name2, i+1)
                #NOTE: this is actually a reduction problem?
                #NOTE: used to be the below code without all these ifs.
                '''df.loc[i1, i2] = res
                res =  (d2[name1][i] * d2[name2][i+1]).sum() /  M[i]'''
                if name2 == 'M':
                    _nexts = list(dropwhile(lambda x: not np.any(x), d2[name2][i+1:]))
                    if not _nexts:
                        res = 0
                    else:
                        #res =  (d2[name1][i] * d2[name2][i+1]).sum() /  M[i]
                        if name2 == 'M' and max(df.loc[i1, 'I%s'%(i+1)], df.loc[i1, 'D%s'%(i+1)]) > 0:
                            res = 1 - max(df.loc[i1, 'I%s'%(i+1)], df.loc[i1, 'D%s'%(i+1)])
                        else:
                            res =  (d2[name1][i] * _nexts[0]).sum() /  M[i]
                            idx = (i+1) + (len(d2[name2][i+1:]) - len(_nexts))
                            i2 = '%s%s' % (name2, idx)
                else:
                    res =  (d2[name1][i] * d2[name2][i+1]).sum() /  M[i]
                print name1, name2, res,    i1, i2
                #if df.ix[i1][i2] ==0:
                if name2 not in selected:
                    df.loc[i1, i2] = res
                    selected.add(name2)
                #NOTE: can't asign using double slice!
                #if df.ix[i1][i2] ==0: df[i1][i2] = res
                #df['%s%s' %(name1, i)]['%s%s' % (name2, i+1)] = res
    for i in range(is_insert.sum(), N, -1):#[::-1]:
        df = df.drop('M%s'%i).drop('D%s'%i)
    def eval(x):
        return 0 if not x[1:].isdigit() else int(x[1:])
#    swith = lambda x: lambda y: y.startswith(x)
    def get(ltr):
        return max( filter(_[0] == ltr, df.columns), key=eval)
    ix, dx, mx = get('I'), get('D'), get('M')
    df.loc[ix, 'end'] = df.loc[dx, 'end'] = df.loc[mx, 'end'] = 1.0
    #TODO: create start probabilities
    return df


makeC = F(np.array) << F(map, np.array) << F(map, list) << F(filter, str.strip) << F(map, str.strip) << str.split
def add_denom(frac):
    return Ratio(frac.numerator, frac.denominator + 1)
def hmmprofile(C, theta):
    N = C.shape[0]
    seqlen = C.shape[1]
    df = makedf(seqlen+1) #, Ratio(0, 0))
    df[:] = Ratio(0, 0)
    is_insert = ((C == EMPTY).sum(axis=0) / N ) > theta
    def update(p, s):
        current = df.loc[p, s]
        df.loc[p, s] = Ratio(current.numerator + 1, current.denominator)
        df.loc[p] = map(add_denom, df.loc[p])
    prev = 'S'
    for seq in C:
        for i, ltr in enumerate(seq):
            off = is_insert[i:].sum()
            if is_insert[i]:
                if ltr != '-': state = st('I', i)
                else: continue
            if not is_insert[i]:
                state = st('D', i+1) if ltr == '-' else st('M', i+1)
            update(prev, state)
            prev = state
        update(prev, 'E')
        prev = 'S'
        l = {'M' : 0, 'D' : 1, 'I' : 2, 'S' : -1, 'E' : 9999999}.__getitem__
    def eval(x):
        #return -1 if not x[1:].isdigit() else (10 * int(x[1:])) + l(x[0])
        return l(x) if len(x) ==1  else (10 * int(x[1:])) + l(x[0])


    cols = sorted(df.columns, key=eval)
    #df = df[cols]
    return df[cols].applymap(float)




def hmmprofile(C, theta):
    N = C.shape[0]
    seqlen = C.shape[1]
    alphabet = set( C.ravel() ) - set(['-'])
    #df = makedf(seqlen+1) #, Ratio(0, 0))
    is_insert = ((C == EMPTY).sum(axis=0) / N ) > theta
    df = makedf((~is_insert).sum()+1) #, Ratio(0, 0))
    EM = pd.DataFrame(0, columns=sorted(alphabet), index=df.index)
    df[:] = EM[:] = Ratio(0, 0)
    def update(df, p, s):
        current = df.loc[p, s]
        df.loc[p, s] = Ratio(current.numerator + 1, current.denominator)
        df.loc[p] = map(add_denom, df.loc[p])
    def update_em(df, state, ltr):
        if ltr != '-':
            df.loc[state] = map(add_denom, df.loc[state])
            current = df.loc[state, ltr]
            df.loc[state, ltr] = Ratio(current.numerator + 1, current.denominator)

    prev = 'S'
    for seq in C:
        for i, ltr in enumerate(seq):
            off = is_insert[:i].sum()
            if is_insert[i]:
                if ltr != '-': state = st('I', i - off)
                else: continue
            if not is_insert[i]:
                state = st('D', i+1 - off) if ltr == '-' else st('M', i+1 - off)
            update(df, prev, state)
            update_em(EM, state, ltr)
            prev = state
        update(df, prev, 'E')
        #update_em(EM, 'E', ltr)
        prev = 'S'
        l = {'M' : 0, 'D' : 1, 'I' : 2, 'S' : -1, 'E' : 9999999}.__getitem__
    def eval(x):
        #return -1 if not x[1:].isdigit() else (10 * int(x[1:])) + l(x[0])
        return l(x) if len(x) ==1  else (10 * int(x[1:])) + l(x[0])


    cols = sorted(df.columns, key=eval)
    #df = df[cols]
    df= df[cols].applymap(float)
    idx = sorted(df.index, key=eval)
    return df.loc[idx], EM.loc[idx].applymap(float)


raw = '''
CBA
CB-
C--
CCA'''
C =  makeC(raw)
df1 = hmmprofile(C, theta=0.289)


raw = '''
E-D-ACA-B
EE--AAAC-
ACD--CBC-
E--B---B-
EAC-AC-CB
DBDDACACB
-CDDACA--
ECD--CA--'''
C =  makeC(raw)
df2 = hmmprofile(C, theta=0.28)

#raw = '''
#DCDABACED
#DCCA--CA-
#DCDAB-CA-
#BCDA---A-
#BC-ABE-AE
#'''
#C =  makeC(raw)
#df3 = hmmprofile(C, theta=0.28)
#raw = '''\nEBA\n E-D\n EB-\n EED\n EBD\n EBE\n E-D\n E-D\n '''
#C =  makeC(raw)
#theta=.289
#df = make_hmm(C, theta=0.289)
#df = hmmprofile(C, theta=0.289)

#raw='''
#ECAEDCEED
#A-AEDC--A
#ECA-DEA-A
#ECAED--EA
#ECACDCA-A
#EDAEDC-EC
#ECAEABAEA
#BBAEDDA--'''
#theta=0.272
raw='''
DCDABACED
DCCA--CA-
DCDAB-CA-
BCDA---A-
BC-ABE-AE
'''
theta=.252
#theta=0.363
raw='''
D-ECB-ABE
D--CBAACE
D--CB-DAE
-EECBAACE
DBECBAAC-
DECBBAECD
DEBEBAEED'''

theta=.333
C =  makeC(raw)
opts=[
('display.float_format',lambda x: ('%.3f' % x) if x != 0 else str(0),),
('display.max_columns', None,),
('display.max_info_rows', 1690785,),
('display.expand_frame_repr', False,),
('display.chop_threshold', None,),]
for a, b in opts:
     pd.set_option(a, b)
#starmap(pd.set_option, opts)






import re
df, em =  hmmprofile(C, theta=theta)
res= '\n--------\n'.join([str(df), str(em)])
res = re.sub(' +', '\t', res)
#with open('e', 'w') as out: out.write(res)
assert res == open('e').read()
print 'passed'

raw='''
A-A
ADA
ACA
A-C
-EA
D-A'''
theta=.358
df, em =  hmmprofile(C, theta=theta)
res= '\n--------\n'.join([str(df), str(em)])
res = re.sub(' +', '\t', res)
print res


#NOTE: pseudocounts
'''
(s + 0.01) / (s + 0.01).sum()
array([ 0.8184466 ,  0.17184466,  0.00970874])
s = df.loc['S'][1:4]
'''
