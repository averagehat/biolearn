


res = np.array([[0.992, 0.988, 0.5, 0.012, 0.008], [0.008, 0.012, 0.5, 0.988, 0.992]])
data = np.array([-3, -2, 0, 2, 3])
centers = np.array([ -2.5, 2.5])

make_dists = partial(dist_matrix, centers.T)
dists = np.apply_along_axis(make_dists, 1, A)


responsibility(dist_matrix(centers, data.T))
