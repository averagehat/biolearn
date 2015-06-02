from func import *
from fn.iters import accumulate
from fn import _, F
from fn.func import flip
import operator as op
from fractions import Fraction
'''
GC-content as an instance of accumulate
end-trimming using dropwhile
website to host bioinformatics codesharing
'''
def hamming(s1, s2):
    assert len(s1) == len(s2)
    return ilen(ifilter(bool, imap(op.ne, s1, s2)))



def gc(A):
    (G+C)/(A+T+G+C) * 100
#gc = parital2(dict(G=1,C=1).get, 0)
gc = F(flip(dict(G=1,C=1).get), 0)
#content = starcompose(op.add, cmp2(gc, gc))
content = compose(op.add, gc)
#content = gc(_) + gc(_)
def content(acc, x):
    return acc + gc(x)
s = 'AGCAACGT'
def content(acc, x):
    return Fraction(acc.numerator + x, acc.denominator+1)
#NOTE: this version of accumulate won't let you pass an initial value.
print tuple(accumulate(imap(gc, s), content))





