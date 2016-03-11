from func import *
import matplotlib
import operator
import numpy as np
#from fn.iters import accumulate
from fn import _, F
X = _
from fn.func import flip
import operator as op
from fractions import Fraction
from matplotlib import pyplot as plt
'''
GC-content as an instance of accumulate
end-trimming using dropwhile
website to host bioinformatics codesharing
'''

class Ratio(Fraction):
    ''' new is for immutable classes.'''
    def __new__(cls, num, denom):
        '''cls is like self really '''
        if denom ==0 and num == 0:
            self = super(Fraction, cls).__new__(cls)
        else:
            self = super(Fraction, cls).__new__(cls)
        self._numerator, self._denominator = num, denom
        return self
    #TODO: fix this for multiplication, etc which fails on denom=0 also unbreak less-than, etc.
    def __float__(self):
        if self._denominator == 0:
            return 0.0
        else:
            return super(Fraction, self).__float__()

Functional programming principles can be productively applied even in small cases. The predictability and compositional capacity of small functions 
allows for modular, extensible code. First-order functions can be used in a variety of generic higher-order functions or "combinators". 
In the future, we will look at how static typing combined with immutability can make this code much safer as an additional benefit.

Let's use calculating [GC-content](https://en.wikipedia.org/wiki/GC-content) as an example.

Where `seq` is some DNA sequence--just a string of characters composed soley of items in the set ```{'A', 'C', 'G', 'T'}```
```python
def gc(seq):
    '''(G+C)/(A+T+G+C) * 100'''
    return seq.count('G') + seq.count('C') / float( len(seq) )
```
from collections import Counter
def gc_counter(seq):
    charmap = Counter(seq)
    #counter is a defaultdict, so this will work even if 'G' or 'C' are not present
    return (charmap['G'] + charmap['C'])/float(sum(charmap.values()))

'''This is all well and good, but what if we want to know the GC-content at a certain place in the read?'''
seq = 'AGCTTAGGCCTTTAAAACCGGGGCCCCCGGAAGCGACTT'
print gc_counter(seq[:10])
''' That works, but will get tiresome and inefficient...'''
gc_counter(seq[0]), gc_counter(seq[:10]), gc_counter(seq[-1])
''' What if we want to create a histogram of the gc-content at each position in the read?'''
N = len(seq)

[gc_counter(seq[:i]) for i  in xrange(1, N)]
''' we have to re-calculate the GC-ratio each time! In terms of runtime, this is equivalent to: '''
sum(i for i in xrange(N) )
''' what happens if we have a longer sequence?'''
sum(i for i in xrange(100) )
sum(i for i in xrange(200) )
sum(i for i in xrange(1000) )
''' Youch! That looks quadritic to me. Method is O(N^2) (about (N^2)/2, where N=length of sequence).'''
gc = F(flip(dict(G=1,C=1).get), 0)
#content = starcompose(op.add, cmp2(gc, gc))
#content = compose(op.add, gc)

#@logcall
#    #NOTE: the problem is that the Fraction class runs gcd to reduce the fraction,making it impossible to represent 0/X where X>1.
#    #additionally, it is impossible to create a fraction 0/0, which would be our starting acc.

#NOTE: the itertools version of accumulate won't let you pass an initial value.

'''
use accumulate for a histogram
use partition to histogram with colors
'''
l = [0, 2, 40, 90, 30, 92, 37, 30]
#def accumulate(iterable, func=operator.add, start=None):
#    it = iter(iterable)
#    if start is None: total = next(it)
#    else: total = start
#    #could also skip yielding total
#    yield total
#    for element in it:
#        total = func(total, element)
#        yield total

def accumulate(iterable, func=operator.add, start=None):
    '''re-implementation of accumulate that allows you to specify a start value
    like with reduce.'''
    'Return running totals'
    # accumulate([1,2,3,4,5]) --> 1 3 6 10 15
    # accumulate([1,2,3,4,5], operator.mul) --> 1 2 6 24 120
    it = iter(iterable)
    if start is None:
        total = next(it)
        yield total
    else:
        total = start
    #could also skip yielding total
    for element in it:
        total = func(total, element)
        yield total

    #if type(acc) == str: return Ratio(val + gc(acc), 2)
def content(ratio, nt):
    '''the advantage here is a rolling, always-true representation of GC-content.'''
    val = 1 if nt in 'GC' else 0
    return Ratio(ratio.numerator + val, ratio.denominator+1)

'''using our function `content` with reduce will give us the total GC-ratio for a given sequence.
`reduce` emulates recursion by computing a new value for each element in the list and passing it onto
the next call with the element as a paramter. '''
#NOTE: don't ever use this code (especiall `_reverse`, which copies and re-creates a new list every "fold"!
def _sum(acc, elem):
    return elem + acc
''' under the hood, this might look like:
[1, 2, 3], 0
-> [1, 2], 3
-> [1],   5
>>>6
'''
def _reverse(acc, elem):
    return [elem] + acc
reduce(_reverse, l, [])
'''
[1, 2, 3, 4], []
-> [1, 2, 3], [4]
-> [1, 2],    [4, 3]
-> [1],       [4, 3, 2]
>>>[4, 3, 2, 1]
'''
'''
notice that at any point during the traversal, the accumulated value
is correct for the traversed part of the list. We'll use that to our advantage later.
Using this model of folding an accumulating paramter over a sequence,
we can model GC-content as a "rolling ratio" over a given sequence of nucleotides.
At any point during the traversal, The ratio will be correct, and the result of "reducing"
the sequence with this model will give us our total GC-content. The following method is not
the most efficient nor the simplest (it requires building a "Ratio" subclass), but it closely
(and flexibly) models the mathematical formula that defines GC-content in the first place. '''
seq = 'AGCTTAGGCCTTTAAAACCGGGGCCCCCGGAAGCGACTT'
print reduce(content, seq, Ratio(0,0))#, None)
'''
Because our `content` function updates the ratio for each partitioning of the list, we can re-use it
to create histogram! viola:'''
def gc_hist(seq):
    func_result = list(accumulate(seq, content, Ratio(0,0)))
    return func_result
'''or'''
gc_hist = partial(accumulate, func=content, start=Ratio(0,0))
#doplot(result)

#necessary because Ratio class is brokers
asfloats = map(float, gc_hist(seq))
THRESH=.5
#TODO: plot above and below as different colors
above, below = map(list, partition(X<THRESH, asfloats))
above, below = map(list, partition(X[1]<THRESH, enumerate(asfloats)))

#nonzero will return those indices where the element is not zero.
#numpy essentially treats true/false as nonzero.
# this also works with N-dimensional arrays, although directly using the index
# will flatten the array. the first
#take until the GC-content reaches a certain threshold.
#np_takewhile = npres[ (npres < THRESH).nonzero()[0][0] ]
''' let's say we want to get all GC-contents below a certain threshold. '''
below_50 = (X < .50001)
'''equivalent to:'''
lambda_blw_50= lambda x: x < .50001
''' applying this function to each element yields true or false.'''
map(below_50, asfloats)
''' we can use this to give us only those scores below the thershold:'''
filter(below_50, asfloats)
'''one of the advantages of programming in a functional style is that functions are eminently reusable.
Python provides a number of "Functor" functions (functions which take other functions as arguments) to
take advantage of this re-usability. Not only can we filter based on the threshold, we can also groupby,
partition, split, etc. our data based on this (or any other) boolean metric. '''


'''
takewhile returns all elements of the list until the boolean operator fails. so This will give us
'''
QUALTHRESH=38
badqual = X < QUALTHRESH
from itertools import takewhile, dropwhile
#trimend = compose_all(reversed, list, qualtrim, reversed, list)
#trimends = compose(trimend, qualtrim)
#seq[slice(len(seq), 0, -1)]
def trim_read(seqrec, QUALTHRESH=38):
    quality = seqrec.letter_annotations['phred_quality']
    badqual = X[1] < QUALTHRESH
    return dropwhile(badqual, izip(seqrec.seq, quality))
l = [0, 2, 40, 90, 30, 92, 37, 30]
#l[:-len(list(takewhile(badqual, reversed(l))))]
#l[ : ilen(qualtrim(reversed(l)))]
#probably faster with takewhile.
'''http://en.wikipedia.org/wiki/Trimming_%28computer_programming%29#Haskell'''
qualtrim = partial(dropwhile, badqual)
''' time to go crazy! let's try writing python in point-free style--by referring
to input and output in the form of fuction composition, absent of their arguments!'''
trim = compose_all(reversed, list, qualtrim)
trim_ends = compose(trim, trim)
list(trim_ends(l))
''' http://en.wikipedia.org/wiki/Tacit_programming  for more info. '''
print list(takewhile(below_50, asfloats))

''' That's all reasonable, and it's nice to see GC-content representated as a ratio. But is it practical?
Well, it will work, but there is a much more efficient way. Array-wise computations like this--which will become
quite large if we get big reads (or god forbid) a whole contig/genome. Additionally, we may want to scale to viewing
multiple reads at once, ie, as a matrix. We will see how the efficient `numpy` library can achieve similar results,
and how we can use these same functional principals--recursion, reduction, and accumulation--to get more leverage (and cleaner code)
out of numpy.'''
npseq =np.array(list(seq))
'''np.sum is actually a specific (read: partial application) of np.reduce!
https://github.com/numpy/numpy/blob/a9c810dd1d8fc1e3c6d0f0ca6310f41795545ec9/numpy/core/_methods.py '''
gccon = ((npseq == 'C') | (npseq == 'G')).sum()/float(len(npseq))
npresult = ((npseq == 'C') | (npseq == 'G')).cumsum()
''' cumulative sum rolls like accumulate. '''
''' divide by the index (starting at one) to simulate the ratio.'''
gcs = ((npseq == 'C') | (npseq == 'G'))
#equivalent to func_result
npres = gcs.cumsum()/np.arange(1,len(npresult)+1, dtype=float)
np_filter_idx = (npres >= .5).nonzero()
'''see np.ufunc.reduce'''
def np_gc_hist(seq):
     npseq = np.array(list(seq))
     gcs = ((npseq == 'C') | (npseq == 'G'))
     npres = gcs.cumsum()/np.arange(1,len(gcs)+1, dtype=float)
     #npres = gcs.cumsum()/xrange(1.0,len(npresult)+1) #, dtype=float)
     return npres


''' using an operator directly (like <) is probably better than a vectorized funciton,
because numpy will handle it internally with C code. A simple work-around:'''
from functools import partial
below_50 = partial(operator.gt, .50001)
'''The `partial` function, similar to the lambda expression above, returns a new function.
The new function is *partially applied*--meaning that we are creating a new function with the
first argument predefined (in this case, .50001). This function takes one less argument. This
works similary to a function defined in a closure, and is equivalent to: '''
def belownum(num):
    def lessthan(other):
        return num > other
    return lessthan
''' In this case, the original `num` argument is available to the returned `lessthan` function
within the closure scope. The result is a curry-able function.'''
belownum(10)(9)
belownum(0)(70)
partial(operator.gt, 10)(9)

print list(takewhile(below_50, asfloats))
def numpy_takewhile(A, boolfunc):
    ''' return array A until the first instance where boolfunc is false.'''
    failing_idxs = np.argwhere(~boolfunc(A))
    if len(failing_idxs) == len(A):
        ''' there is no matching element in the array. '''
        return np.array([])
    ''' return array A until the first istnace (index) of failure.'''
    return A[: failing_idxs[0]]

print numpy_takewhile(npres, below_50)
def doplot(result):
    fig = plt.figure()
    fig.set_size_inches( 20.0, 8.0 )
    gs = matplotlib.gridspec.GridSpec(1,2, width_ratios=[20,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax1.plot(result)
    '''draw a line for the threshold.'''
    ax1.axhline(y=THRESH, color='green')
    fig.show()

#fig.savefig('red')
#ax.plot(*zip(*above), color='b')
#ax.plot(*zip(*below), color='r')
''' objects are badkay '''
def zipqual(seqrec):
    return izip(seqrec.seq, seqrec.letter_annotations['phred_quality'])
_badqual = compose(badqual, X[1])
_qualtrim = partial(dropwhile, _badqual)
_trim = compose_all(reversed, list, _qualtrim)
_first_trim = compose(_trim, zipqual)
_trim_ends = compose(_trim, _first_trim)


f = '/home/AMED/michael.panciera/projects/data/1.R1.unmap.fastq'
from Bio import SeqIO
recs = SeqIO.parse(open(f), 'fastq')
''' returns an iterator. '''

''' timing stuff'''
from random import choice
N = 1000000
alpha = list('AGCT')
seq = ''.join( choice(alpha) for _ in xrange(N) )
'''or '''
from fn.iters import repeatfunc
seq = ''.join(repeatfunc(F(choice, alpha), N))
gcl = compose(list, gc_hist)
'''
In [257]: %timeit res = np_gc_hist(seq)
1000 loops, best of 3: 1.05 ms per loop
'''

'''
#NOTE: without actually evaluating anything! (pure iterators)
In [258]: %timeit funcres = gc_hist(seq)
1000000 loops, best of 3: 417 ns per loop
'''
'''

In [265]: %timeit funcres = list(gc_hist(seq))
1 loops, best of 3: 971 ms per loop


In [267]: %timeit funcres = gcl(seq)
1 loops, best of 3: 987 ms per loop
'''
