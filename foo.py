from itertools import repeat
from func import apply_each, compose
from functools import partial
from itertools import repeat
import operator as op
from utils import slider
import numpy as np
#compute distance for a ll centers from all points
#transofrm matrix via
#compute further by setting val/sum(vals)

#val = 1/(distnace)**2

def rowdiv(row):
    return row/row.sum()


def make_center(row, ones, data):
    return np.dot(row, data)/np.dot(row, ones)



make_np = compose(np.array, list)


norm_matrix = partial(np.linalg.norm, axis=1)
make_dists = compose(norm_matrix, op.sub)
centers = np.array([[-2.5], [2.5]])
nctrs = len(centers)
data = np.array([-3, -2, 0, 2, 3])
#hm = np.empty( (len(centers), len(data)) )
def p_part(A, B, stiffness=0.5):
    if B.shape[0] == 1:
        dist = op.sub
    return np.e**(-stiffness * dist(A,B))


def get_hm(data, centers):
    data_grid = np.array(list(repeat(data, nctrs))).T
    intermed =  1/(data_grid-centers.T)**2
    #intermed =  p_part(data_grid, centers.T)
    #intermed = np.e**(-.5* (data_grid-centers.T))
    hm = np.apply_along_axis(rowdiv, 1, intermed).T
    return hm


def get_centers(hm, data):
    #n_centers = np.empty( (k, data.shape[1]) )
    ones = np.ones(hm.shape[1])
    row_funcs = [partial(make_center, row, ones) for row in hm]
    _all = [partial(np.apply_along_axis, row_func, 0) for row_func in row_funcs]
    #for row in hm:
    res= make_np(apply_each(_all, data))
    return res
#np.apply_along_axis(row_funcs, 0, data)
#new_centers = make_np( apply_each(row_funcs, data) )



hm = get_hm(data, centers)
print hm

#n_centers = get_centers(hm, data)
#print n_centers






from itertools import product
from fn import F
import operator as op
s = 'ACGCGGCTCTGAAA'
k = 2
kmers = sorted(map(''.join, product(*([list('ACGT')]*k))))
cnt_with_overlap = lambda x: sum(map(lambda y: x==y, slider(s, len(x))))
print map(cnt_with_overlap, kmers)

sym_num = {'A':0,'C':1,'G':2,'T':3}.__getitem__
def sym(acc, c):
    return sym_num(x) + sym_num(y)

sum(map(sym_num, 'AGT'))



def ptn(s):
    if s == '': return 0
    else:
        return (4*ptn(s[:-1])) + sym_num(s[-1])
