import unittest
from numpy.testing import assert_array_almost_equal as aa
import numpy as np
from numpy import nan
from npgraph import djs, nbr_idx, uniq_path
A = np.array
class TestGraph(unittest.TestCase):


    def setUp(self):
        pass

    def test_aanbr_idx(self):
        _In = np.array([[0, nan, 2],
                        [nan, 0, 5],
                        [2,   5, 0]])

        aa(nbr_idx(_In, 0), A([2]))
        aa(nbr_idx(_In, 1), A([2]))
        aa(nbr_idx(_In, 2), A([0, 1]))


    def test_alg_n1(self):
        _In = np.array([[0, 2],
                        [2, 0]])
        aa(djs(_In)[0], A([nan, 0]))
        aa(djs(_In)[1], A([0, 2]))
    
    def test_alg_n2(self):
        _In = np.array([[0, nan, 2],
                        [nan, 0, 5],
                        [2,   5, 0]])

        aa(djs(_In)[0], A([nan, 2, 0]))
        aa(djs(_In)[1], A([0, 7, 2]))


class TestUniqPath(unittest.TestCase):
    def test_alg_n2_uniq(self):
        _In = np.array([[0, nan, 2],
                        [nan, 0, 5],
                        [2,   5, 0]])
        actual = uniq_path(_In, 2)
        self.assertListEqual(actual, [0, 2] )
        actual = uniq_path(_In, 1)
        self.assertListEqual(actual, [0, 2, 1] )

if __name__ == '__main__':
    unittest.main()


