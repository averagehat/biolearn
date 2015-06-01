from func import compose, partial2
import numpy as np
from functools import partial
import operator as op
from itertools import islice, ifilter
from fn.iters import splitby
from operator import methodcaller
from func import compose_all, pmap, psplit, _not


norm_matrix = partial(np.linalg.norm, axis=1)
dist_matrix = compose(norm_matrix, op.sub)
min_dist_position = compose(np.argmin, dist_matrix)

def gravity(M):
    return M.sum(axis=0)/float(len(M))

def k_means_cluster(data, centers=None, k=None):
    assert (centers is not None or k is not None)
    if centers is None:
        centers = np.empty((k, data.shape[1]))
        centers[:] = data[:k]
    old_centers = centers.copy()

    mincenters = partial(min_dist_position, centers)
    data_by_cluster = np.apply_along_axis(mincenters, 1, data)
    for i, _ in enumerate(centers):
        cluster = data[data_by_cluster == i]
        centers[i] = gravity(cluster)
    if np.allclose(old_centers, centers):
        return centers, data_by_cluster
    else:
        return k_means_cluster(data, centers) #, runs=runs+1)

stiffen = lambda x: 1/x**2
row_sum = partial(np.sum, axis=0)
col_sum = partial(np.sum, axis=1)



def responsibility(distances):
    ''' return a KxN matrix'''
    stiffs = stiffen(distances)
    return stiffs.T/stiffs.sum(axis=1)



def make_centers(res, data):
    onematrix = np.ones(res.shape[1])
    return np.dot(res, data).T/np.dot(res, onematrix)

def soft_k_means_cluster(data, centers=None, k=None, step=0):
    if centers is None:
        centers = np.empty((k, data.shape[1]))
        centers[:] = data[:k]
    make_dists = partial(dist_matrix, centers)
    dists = np.apply_along_axis(make_dists, 1, data)
    #KxN matrix
    respons = responsibility(dists)
    print respons
    centers = make_centers(respons, data)
    print centers

    #data_by_cluster = np.apply_along_axis(mincenters, 1, data)
    for i, _ in enumerate(centers):
        cluster = data[data_by_cluster == i]
        centers[i] = gravity(cluster)
    if steps==100:
        return centers, data_by_cluster
    else:
        return k_means_cluster(data, centers, step=step+1) #, runs=runs+1)




#def distortion(data, centers):
#    mindests = np.array(map(partial(mindist, centers=centers), data))
#    return (mindests**2).mean()
#
#def mindist_coordinates(d, centers):
#    ddist = partial(dist, d)
#    return ddist(centers).argmin()
#    #return min(centers, key=ddist)


def makematrices(s):
    _centers, _data = splitby(_not(isin('------')), ifilter(bool, s))
    #centers = map(makenp, islice(_centers, 1, None))
    #data = map(makenp, islice(_data, 1, None))
    centers = makenp(islice(_centers, 1, None))
    data = makenp(islice(_data, 1, None))
    return centers, data



isin = partial(methodcaller, '__contains__')
makearray = compose_all(np.array, pmap(np.array), pmap(float), psplit(' '))
makenp = compose(np.array, pmap(makearray))
def get_in_out(s):
    raw_in, raw_out = splitby(_not(isin('Output')), ifilter(bool, s))
    k = int(next(raw_in).split(' ')[0])
    _in = makenp(raw_in)
    _out =makenp(islice(raw_out, 1, None))
    return _in, _out, k


lines = open('Lloyd.txt').readlines()
input, expected, k = get_in_out(lines)
print soft_k_means_cluster(input, k=3)


from matplotlib import pyplot
colors = np.random.rand(k)
centers, clusters = k_means_cluster(input[:,:2], k=k)

pyplot.scatter(input[:,0], input[:,1], c=clusters)
pyplot.scatter(centers[:,0], centers[:,1], marker='x', s=300, linewidth=4, c='black')
pyplot.show()



