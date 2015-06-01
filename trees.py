from numpy import nan
import numpy as np
from operator import itemgetter
import re
from numpy.testing import assert_almost_equal
import unittest
from itertools import ifilter, product, ifilterfalse
from func import compose, typecheck, ilen, starcompose, reverse
from numbers import Number as num


filterfst = compose(next, ifilter)

class TestTrees(unittest.TestCase):

  def setUp(self):
    self.complex='''
    0->4:11
    1->4:2
    2->5:6
    3->5:7
    4->0:11
    4->1:2
    4->5:4
    5->4:4
    5->3:7
    5->2:6
    '''
    self.M, self.V = getdists(self.complex)


  def test_alg(self):
        M = np.array([ [0, 3, 5  ],
                       [3, 0, nan],
                       [5, nan, 0] ])
        actual = mindist(M, np.full( (3, 3), False), 1, 2)[-1]
        self.assertEquals(actual, 8)
  def test_real_data_leaves(self):
    actual = get_leaves(self.M)
    self.assertListEqual(actual, [0, 1, 2, 3])

  def test_real_data_result(self):
      actual = mindist(self.M, self.V, 0, 1)[-1]
      self.assertEquals(actual, 13)



  def test_simple(self):
      _in='''
      0->1:3
      1->0:3
      '''
      expected  = np.array( [ [0, 3],
                              [3, 0] ] )
      actual = getdists(_in)[0]
      assert_almost_equal(actual, expected)

  def test_get_leaves(self):
        M = np.array([ [0, 3, 5  ],
                       [3, 0, nan],
                       [5, nan, 0] ])
        actual = get_leaves(M)
        self.assertListEqual(actual, [1, 2])

  def test_get_leaves_empty(self):
        M2  = np.array( [ [0, 3, 1], [3, 0, 1], [2, 0, 2 ] ])
        self.assertListEqual(get_leaves(M2), [])


  def test_with_three_nodes(self):
      _in='''
      0->1:3
      1->0:3
      2->0:5
      0->2:5
      '''
      expected = np.array([ [0, 3, 5  ],
                            [3, 0, nan],
                            [5, nan, 0] ])
      actual = getdists(_in)[0]
      assert_almost_equal(actual, expected)

  def test_full_matrix_compute(self):
    print self.M
    expected=np.array([[  0 ,   13,  21,  22],
                       [  13,   0 ,  12,  13],
                       [  21,   12,  0 ,  13],
                       [  22,   13,  13,  0 ] ])

    actual = compute(self.complex)
    assert_almost_equal(actual, expected)

def t(tup):
    return tup[-1]



@typecheck(np.ndarray, np.ndarray, num, num)
def mindist(M, v, start, end, row_wise=False, steps=0):

    i, j = start, end

    if M[i, j] == 0: print 'zerohero';return None
    if steps > 1000: import sys; sys.exit(0)
    if not np.isnan(M[i, j]):
        return (i, j, M[i, j])
    if v[start, end]: return None #print 'visited'; return None
    v[start, end] = v[end, start] = True
    getpath = lambda x, M=M, end=end: mindist(M, v, end, x, row_wise=(not row_wise), steps=steps+1)
    row_neighbors = range(0, start) + range(start+1, M.shape[0])
    col_neighbors = range(0, end) + range(end+1, M.shape[1])
    row_neighbors.remove(end)
    col_neighbors.remove(start)
    neighbors = row_neighbors if row_wise else col_neighbors
    res = filter(bool, map(getpath, neighbors))
    nexti, nextj, next_dist = min(res, key=t)
    M[nexti, nextj] = M[nextj, nexti] = next_dist
    part =  mindist(M, v, start, nextj, row_wise=row_wise, steps=steps+1)
    if not part: return None
    M[start, nextj] = M[nextj, start] = part[2]
    print part
    return part[0], part[1], part[2] + next_dist






def get_leaves(M):
    return [i for i, row in enumerate(M) if ilen(ifilterfalse(np.isnan, row)) == 2]

def parse_line(line):
    '''
    :return (fromNode, toNode, distance)
    '''
    reg = re.compile('\d+')
    fromNode, toNode, distance = map(int, reg.findall(line))
    return fromNode, toNode, distance

def getdists(raw):
    ''' return 2d matrix '''
    lines = ifilter(str.strip, raw.splitlines())
    codes = map(parse_line, lines)
    n = len(set(map(itemgetter(0), codes)))
    dists = np.full( (n, n), np.nan)
    visits = np.full( (n, n),False)
    for fromNode, toNode, distance in codes:
        dists[fromNode, toNode] = distance
    np.fill_diagonal(dists, 0)
    return dists, visits

def leave_dists(M, V):
    leaves = get_leaves(M)
    ldists = np.full( (len(leaves), len(leaves) ) , 0)
    c_prod = product(leaves, leaves)
    uniques =  uniq_prods(c_prod)
    for cell in uniques:
        dist = mindist(M, V, cell[0], cell[1])[-1]
        V =  np.full( V.shape,False)
        ldists[cell] = ldists[reverse(cell)] = dist
    return ldists

compute = starcompose(leave_dists, getdists)

def uniq_ignore_order(acc, x):
     if reverse(x) in acc or (len(x) > len(set(x))):
         return acc
     return set([x]) | acc

uniq_prods = lambda seq: reduce(uniq_ignore_order, seq, set())


def get_products(ints, times=2):
    return uniq_prods(product(*([ints]*times)))




