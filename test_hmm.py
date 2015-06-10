from hmm import viterbi, backtrk, od, prfull, forward, _toints, pr_hidden_path
import numpy as np

a = np.array([[ 0.377,   0.623], [0.26,    0.74]])
p = 'ABABBBAAAA'


assert (0.000384928691755 -   pr_hidden_path(a, p))  < 0.000001


p = 'BBABBBABBAABABABBBAABBBBAAABABABAAAABBBBBAABBABABB'
a= np.array([[     0.863 ,  0.137   ], [     0.511 ,  0.489]])

assert (3.26233331904e-21 -   pr_hidden_path(a, p)  ) < 0.000001
p = 'BBABAABABABABABBAAABABBBAAAAAABBBBBBBABABBAABAABAB'

a = np.array([[0.339,    0.661   ], [  0.64,    0.36]])

AE = np.testing.assert_almost_equal
AE(prfull('datas/foo'), 3.42316482177e-35)
TM = np.array([[0.641,     0.359 ],
               [0.729,     0.271 ]])
EM = np.array([[0.117,     0.691,   0.192],
               [0.097,     0.42 ,   0.483]])
ems = 'xyxzzxyxyy'
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
A = np.array
s = 'GGCA'
tm = A( [ [0.5, .5],
          [.4, .6]])
em = A([[ 0.2, .3, .3, .2],
       [ 0.3, .2, .2, .3]])
score = forward(s , tm, em, _toints)
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
