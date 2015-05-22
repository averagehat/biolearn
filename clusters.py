import numpy as np
from func import *
from math import sqrt
#Center = namedtuple("Center", ["point", "matrix"])

def dist(c, d):
    return sqrt(np.square(c-d).sum())
#@typecheck(np.matrix, dict(centers=tuple))
def get_max_idx_excluding_other(matrix, dists, other):
    _max = dists.argmax()

    if tuple(matrix[_max]) in map(tuple, other):
        dists.mask[_max] = True
        return get_max_idx_excluding_other(matrix, dists, other)
    return _max



def compute(matrix, centers, p, k, i=0):#, calls=0):
    #
    assert i <= k
    if len(centers) == k:
        return centers
    #dist_makers = map(partial(partial, dist), matrix)
    #raw_dists = apply_each(dist_makers,  centers)
    #raw_dists = map(partial(apply_each, dist_makers),  centers)
    #print list(raw_dists)
    #distances = np.ma.array(list(raw_dists), mask=False)
    #for point in matrix:
    dists = np.ma.array([sum(map(partial(dist, point), centers)) for point in matrix], mask=False)
    #dists.data.sort()
    farthest_idx = get_max_idx_excluding_other(matrix, dists, centers)
    p = matrix[farthest_idx]
    return compute(matrix, centers + tuple([p]), p, k, i+1)

    #farthest_idx = dists.argmax()
#    distances = np.ma.array(sum(np.array(map(list, raw_dists))), mask=False)
    #farthest_idx = get_max_idx_excluding_other(matrix, distances, centers)
    #matrix.mask[farthest_idx] = True
    #use a masked array
    #distance_matrix = lambda p: matrix.apply(partial(dist, p))
    #distance_matrix = lambda x, matrix=matrix: partial(dist, x)#(matrix)
    #distance_matrix = lambda x, matrix=matrix: np.vectorize(partial(dist, x))(matrix)
    #return compute(matrix, centers+ (p), p, k, i+1)
#    v_distance_makers = map(partial(partial, v_dist), centers)
#    #v_distance_makers = map(np.vectorize, distance_makers)
#    distances_gen = apply_each(v_distance_makers, matrix)
#    #distnaces = np.fromiter(distances_gen, dtype=float)
#    distances = np.array(list(distances_gen)) #.sum(axis=1)
#    distances = map(partial(np.sum, axis=1), distances)

    #centers[i] = p
from operator import methodcaller
lines = open('FarthestFirstTraversal.txt', 'r').readlines()
k, m = map(int, lines[1].split())
_in, out = partition(methodcaller('__contains__', 'Output'), lines[2:])
matrix = np.ma.array(map(pmap(float), map(str.split, _in)), mask=False)




expected='''
0.8 12.0 17.5 0.9 7.2
0.3 16.4 8.9 34.6 24.6
32.3 1.9 5.1 16.2 8.8
23.1 31.1 3.6 0.8 0.3
'''


'''
return compute(centers | set(p), p, k, calls+1)
there is no numpy apply, only vectorized funcions.
can np.vectorize the function if necessary (for example, it uses "if x")
argmax returns the index of the max value
'''
#centers = np.empty(shape=(k, m))
start = matrix[0]

actual = compute(matrix, tuple([start]), k=k, p=start)#, centers)
print actual
    #farthest = max(sum(starmap(compute_matrix, izip(centers, repeat(p)))), key=sum)

    #distances = centers.map(distance_matrix)
#    farthest = max(sum(starmap(compute_matrix, izip(centers, repeat(p)))), key=sum)
#    return compute(centers | farthest, farthes, k+1)

m=np.array([
[0.0, 0.0],
[5.0, 5.0],
[0.0, 5.0],
[1.0, 1.0],
[2.0, 2.0],
[3.0, 3.0],
[1.0, 2.0]])

start = m[0]

print compute(m, tuple([start]), k=3, p=start)#, centers)














