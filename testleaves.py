import unittest
from numpy.testing import assert_almost_equal
from trees import get_products, getdists
from itertools import starmap, ifilter
from functools import partial
from func import pmap, compose, compose_all
import numpy as np

class TestLimb(unittest.TestCase):

    def test_simple_example(self):
        j = 1
        D=np.array([
        [0 ,  13,  21,  22],
        [13,  0 ,  12,  13],
        [21,  12,  0 ,  13],
        [22,  13,  13,  0]])
        actual = limb_len(D, j)
        self.assertEquals(actual, 2)

    def test_complex_example(self):
        _in = '''0 6806 3415 3666 8467 1175 6105 4705 1537 5183 4463 2616 2156 9275 3315 7970 4217 2632 7561 8857 4047 9129 4972 3729 8378
6806 0 9639 9890 3615 7399 2689 10929 6371 11407 2849 5972 8380 4423 4905 3118 3723 8856 2709 4005 3737 4277 11196 9953 3526
3415 9639 0 1319 11300 2886 8938 2358 4370 2836 7296 5449 2541 12108 6148 10803 7050 1061 10394 11690 6880 11962 2625 1382 11211
3666 9890 1319 0 11551 3137 9189 2035 4621 2513 7547 5700 2792 12359 6399 11054 7301 1312 10645 11941 7131 12213 2302 1059 11462
8467 3615 11300 11551 0 9060 4350 12590 8032 13068 4510 7633 10041 1978 6566 2249 5384 10517 2904 1560 5398 1832 12857 11614 1207
1175 7399 2886 3137 9060 0 6698 4176 2130 4654 5056 3209 1627 9868 3908 8563 4810 2103 8154 9450 4640 9722 4443 3200 8971
6105 2689 8938 9189 4350 6698 0 10228 5670 10706 2148 5271 7679 5158 4204 3853 3022 8155 3444 4740 3036 5012 10495 9252 4261
4705 10929 2358 2035 12590 4176 10228 0 5660 2166 8586 6739 3831 13398 7438 12093 8340 2351 11684 12980 8170 13252 1955 1736 12501
1537 6371 4370 4621 8032 2130 5670 5660 0 6138 4028 2181 3111 8840 2880 7535 3782 3587 7126 8422 3612 8694 5927 4684 7943
5183 11407 2836 2513 13068 4654 10706 2166 6138 0 9064 7217 4309 13876 7916 12571 8818 2829 12162 13458 8648 13730 1033 2214 12979
4463 2849 7296 7547 4510 5056 2148 8586 4028 9064 0 3629 6037 5318 2562 4013 1380 6513 3604 4900 1394 5172 8853 7610 4421
2616 5972 5449 5700 7633 3209 5271 6739 2181 7217 3629 0 4190 8441 2481 7136 3383 4666 6727 8023 3213 8295 7006 5763 7544
2156 8380 2541 2792 10041 1627 7679 3831 3111 4309 6037 4190 0 10849 4889 9544 5791 1758 9135 10431 5621 10703 4098 2855 9952
9275 4423 12108 12359 1978 9868 5158 13398 8840 13876 5318 8441 10849 0 7374 3057 6192 11325 3712 716 6206 1332 13665 12422 2015
3315 4905 6148 6399 6566 3908 4204 7438 2880 7916 2562 2481 4889 7374 0 6069 2316 5365 5660 6956 2146 7228 7705 6462 6477
7970 3118 10803 11054 2249 8563 3853 12093 7535 12571 4013 7136 9544 3057 6069 0 4887 10020 2407 2639 4901 2911 12360 11117 2160
4217 3723 7050 7301 5384 4810 3022 8340 3782 8818 1380 3383 5791 6192 2316 4887 0 6267 4478 5774 1148 6046 8607 7364 5295
2632 8856 1061 1312 10517 2103 8155 2351 3587 2829 6513 4666 1758 11325 5365 10020 6267 0 9611 10907 6097 11179 2618 1375 10428
7561 2709 10394 10645 2904 8154 3444 11684 7126 12162 3604 6727 9135 3712 5660 2407 4478 9611 0 3294 4492 3566 11951 10708 2815
8857 4005 11690 11941 1560 9450 4740 12980 8422 13458 4900 8023 10431 716 6956 2639 5774 10907 3294 0 5788 914 13247 12004 1597
4047 3737 6880 7131 5398 4640 3036 8170 3612 8648 1394 3213 5621 6206 2146 4901 1148 6097 4492 5788 0 6060 8437 7194 5309
9129 4277 11962 12213 1832 9722 5012 13252 8694 13730 5172 8295 10703 1332 7228 2911 6046 11179 3566 914 6060 0 13519 12276 1869
4972 11196 2625 2302 12857 4443 10495 1955 5927 1033 8853 7006 4098 13665 7705 12360 8607 2618 11951 13247 8437 13519 0 2003 12768
3729 9953 1382 1059 11614 3200 9252 1736 4684 2214 7610 5763 2855 12422 6462 11117 7364 1375 10708 12004 7194 12276 2003 0 11525
8378 3526 11211 11462 1207 8971 4261 12501 7943 12979 4421 7544 9952 2015 6477 2160 5295 10428 2815 1597 5309 1869 12768 11525 0
'''
        D = parsematrix(_in)
        actual = limb_len(D, 2)
        self.assertEquals(actual, 534)

def parsematrix(raw):
    _in = filter(bool, raw.split('\n'))
    return  np.ma.array(map(pmap(int), map(str.split, _in)), mask=False)

def parent_len(D, j, i, k):
    return (D[i,j] + D[j,k] - D[i,k])/2

def limb_len(D, j):
    ''' return the length of the limb from some internal node to leaf j.'''
    other_rows = range(0, j) + range(j+1, D.shape[0])
    iks = get_products(other_rows)
    j_parent = partial(parent_len, D, j)
    return min(starmap(j_parent, iks))


from numpy import nan
def test_add_phylo_2D():
      _in = np.array([ [0, 3, 5  ], [3, 0, nan], [5, nan, 0] ])
      expected='''0->1:3
      0->2:5
      '''
      actual = additive_phylo(_in, 2)
      assert expected == actual

filterfst = compose(next, ifilter)
def str_row(D, j):
    row = D[j]
    p = (str(j)+"->{0}:{1}").format
    return '\n'.join(starmap(p, enumerate(row)))

def str_matrix(D):
    d_str = partial(str_row, D)
    return '\n'.join(map(d_str, xrange(D.shape[0])))

# method that gets node with matching distance
def get_match_dst(D, j, dist):
    assert dist != 0
    return (D[j] == dist).argmax()
    #return D[i, (D[i] == dist)]

def get_match_dists(D, j, dist):
    assert dist != 0
    return (D[j] == dist).nonzero()


def non_diag(D, j):
    return  range(0, j) + range(j+1, D.shape[0])

nondiag_products = compose_all(list, get_products, non_diag)
nondiag_products3 = compose_all(list, partial(get_products, times=3), non_diag)
products3 = compose_all(list, partial(get_products, times=3))
def additive_phyloZ(D, n):
    if n == 2:
        return  str_matrix(D)
    ll = limb_len(D, n)
    non_diag = range(0, n) + range(n+1, D.shape[0])
    D[non_diag, j] -= ll
    D[j, non_diag] -= ll
    # get matching i, n, k
    D.mask[n] = D.mask[:, n] = True
    T = additive_phylo(D, n-1)
    v_candidates_i = get_match_dst(T, i, x)
    along_path = lambda c: D[k, c] + D[i, c] == x
    v = filterfst(along_path, v_candidates_i)
    # add leaf n back to T by creating limb (v, n)
    return T


def get_matching_nodes3(D):
    ''' :return (i, k) where n is the center '''
    def has_center_c(ikj, D=D):
        i, k, j  = ikj
        return D[i, k] == D[i, j] + D[j, k]
    candidates = list(get_products(xrange(D.shape[0]), 3))
    return filterfst(has_center_c, candidates)


def get_matching_nodes(D, n):
    ''' :return (i, k) where n is the center '''
    def has_center_c(ik, j=n, D=D):
        i, k = ik
        return D[i, k] == D[i, j] + D[j, k]
    candidates = nondiag_products(D, n)
    return filterfst(has_center_c, candidates)

#E = 2*L - 2
ltrs = ['i', 'j', 'k', 'l']
def additive_phylo(D, n):
    if n <= 1:
        dst = D.data[0,1]
        return {ltrs[0] : [(ltrs[1], dst)], ltrs[1] : [(ltrs[0], dst)] }
    ll = limb_len(D, n)
    nondiag = non_diag(D, n)
    D[nondiag, n] -= ll
    D[n, nondiag] -= ll
    #i, k = get_matching_nodes(D, n)
    i, k, _n = get_matching_nodes3(D)
    #x = D[i, n]
    D.mask[n] = D.mask[:, n] = True
    T = additive_phylo(D, n-1)
    '''add leaf n back to T by creating limb (v, n)'''
    v ='v'+str(D.shape[0] - n)
    #TODO: connect vertices (internal nodes) properly
    #T[v] = [(ltrs[n], ll), (ltrs[i], D.data[i, n]), (ltrs[k], D.data[k, n])]
    T[v] = [(_n, ll), (i, D.data[i, _n]), (k, D.data[k, _n])]
    return T

_in =np.ma.array([[0, 21, 22],
              [21, 0, 13],
              [22, 13, 0]], mask=False)
print additive_phylo(_in, 2)


n=np.ma.array([
[0, 13, 21, 22],
[13, 0, 12, 13],
[21, 12, 0, 13],
[22, 13, 13, 0]], mask=False)
print additive_phylo(n, 3)

D = parsematrix('''0 3036 4777 1541 2766 6656 2401 4119 7488 4929 5344 3516 1485 6392 2066 3216 7008 7206 1187 6491 3379 6262 6153 4927 6670 4997 9010 5793 9032
3036 0 6323 3087 4312 8202 1619 5665 9034 2205 6890 966 3031 7938 3612 492 4284 8752 2023 8037 4925 3538 3429 6473 3946 2273 10556 7339 10578
4777 6323 0 4054 2571 3455 5688 1972 4287 8216 2143 6803 4126 3191 3183 6503 10295 4005 4474 3290 2964 9549 9440 1726 9957 8284 5809 2592 5831
1541 3087 4054 0 2043 5933 2452 3396 6765 4980 4621 3567 890 5669 1343 3267 7059 6483 1238 5768 2656 6313 6204 4204 6721 5048 8287 5070 8309
2766 4312 2571 2043 0 4450 3677 1913 5282 6205 3138 4792 2115 4186 1172 4492 8284 5000 2463 4285 1173 7538 7429 2721 7946 6273 6804 3587 6826
6656 8202 3455 5933 4450 0 7567 3851 1082 10095 1872 8682 6005 1992 5062 8382 12174 800 6353 1047 4843 11428 11319 3053 11836 10163 2604 2545 2626
2401 1619 5688 2452 3677 7567 0 5030 8399 3512 6255 2099 2396 7303 2977 1799 5591 8117 1388 7402 4290 4845 4736 5838 5253 3580 9921 6704 9943
4119 5665 1972 3396 1913 3851 5030 0 4683 7558 2539 6145 3468 3587 2525 5845 9637 4401 3816 3686 2306 8891 8782 2122 9299 7626 6205 2988 6227
7488 9034 4287 6765 5282 1082 8399 4683 0 10927 2704 9514 6837 2824 5894 9214 13006 1172 7185 1879 5675 12260 12151 3885 12668 10995 2144 3377 2166
4929 2205 8216 4980 6205 10095 3512 7558 10927 0 8783 2859 4924 9831 5505 1891 3719 10645 3916 9930 6818 2973 2864 8366 3381 1708 12449 9232 12471
5344 6890 2143 4621 3138 1872 6255 2539 2704 8783 0 7370 4693 1608 3750 7070 10862 2422 5041 1707 3531 10116 10007 1741 10524 8851 4226 1233 4248
3516 966 6803 3567 4792 8682 2099 6145 9514 2859 7370 0 3511 8418 4092 1146 4938 9232 2503 8517 5405 4192 4083 6953 4600 2927 11036 7819 11058
1485 3031 4126 890 2115 6005 2396 3468 6837 4924 4693 3511 0 5741 1415 3211 7003 6555 1182 5840 2728 6257 6148 4276 6665 4992 8359 5142 8381
6392 7938 3191 5669 4186 1992 7303 3587 2824 9831 1608 8418 5741 0 4798 8118 11910 2542 6089 1827 4579 11164 11055 2789 11572 9899 4346 2281 4368
2066 3612 3183 1343 1172 5062 2977 2525 5894 5505 3750 4092 1415 4798 0 3792 7584 5612 1763 4897 1785 6838 6729 3333 7246 5573 7416 4199 7438
3216 492 6503 3267 4492 8382 1799 5845 9214 1891 7070 1146 3211 8118 3792 0 3970 8932 2203 8217 5105 3224 3115 6653 3632 1959 10736 7519 10758
7008 4284 10295 7059 8284 12174 5591 9637 13006 3719 10862 4938 7003 11910 7584 3970 0 12724 5995 12009 8897 1442 2699 10445 1088 3447 14528 11311 14550
7206 8752 4005 6483 5000 800 8117 4401 1172 10645 2422 9232 6555 2542 5612 8932 12724 0 6903 1597 5393 11978 11869 3603 12386 10713 2694 3095 2716
1187 2023 4474 1238 2463 6353 1388 3816 7185 3916 5041 2503 1182 6089 1763 2203 5995 6903 0 6188 3076 5249 5140 4624 5657 3984 8707 5490 8729
6491 8037 3290 5768 4285 1047 7402 3686 1879 9930 1707 8517 5840 1827 4897 8217 12009 1597 6188 0 4678 11263 11154 2888 11671 9998 3401 2380 3423
3379 4925 2964 2656 1173 4843 4290 2306 5675 6818 3531 5405 2728 4579 1785 5105 8897 5393 3076 4678 0 8151 8042 3114 8559 6886 7197 3980 7219
6262 3538 9549 6313 7538 11428 4845 8891 12260 2973 10116 4192 6257 11164 6838 3224 1442 11978 5249 11263 8151 0 1953 9699 1104 2701 13782 10565 13804
6153 3429 9440 6204 7429 11319 4736 8782 12151 2864 10007 4083 6148 11055 6729 3115 2699 11869 5140 11154 8042 1953 0 9590 2361 2592 13673 10456 13695
4927 6473 1726 4204 2721 3053 5838 2122 3885 8366 1741 6953 4276 2789 3333 6653 10445 3603 4624 2888 3114 9699 9590 0 10107 8434 5407 2190 5429
6670 3946 9957 6721 7946 11836 5253 9299 12668 3381 10524 4600 6665 11572 7246 3632 1088 12386 5657 11671 8559 1104 2361 10107 0 3109 14190 10973 14212
4997 2273 8284 5048 6273 10163 3580 7626 10995 1708 8851 2927 4992 9899 5573 1959 3447 10713 3984 9998 6886 2701 2592 8434 3109 0 12517 9300 12539
9010 10556 5809 8287 6804 2604 9921 6205 2144 12449 4226 11036 8359 4346 7416 10736 14528 2694 8707 3401 7197 13782 13673 5407 14190 12517 0 4899 1758
5793 7339 2592 5070 3587 2545 6704 2988 3377 9232 1233 7819 5142 2281 4199 7519 11311 3095 5490 2380 3980 10565 10456 2190 10973 9300 4899 0 4921
9032 10578 5831 8309 6826 2626 9943 6227 2166 12471 4248 11058 8381 4368 7438 10758 14550 2716 8729 3423 7219 13804 13695 5429 14212 12539 1758 4921 0''')
print additive_phylo(D, 28)

D2 = parsematrix('''0   295 306 497 1081    1091    1003    956 954
295 0   309 500 1084    1094    1006    959 957
   306 309     0   489 1073    1083    995 948 946
   497 500 489 0   1092    1102    1014    967 965
1081    1084    1073    1092    0   818 1056    1053    1051
1091    1094    1083    1102    818 0   1066    1063    1061
  1003    1006    995 1014    1056    1066    0   975 973
  956 959 948 967 1053    1063    975 0   16
  954 957 946   965 1051    1061    973 16  0''')
#print additive_phylo(D2, D2.shape[0]-1)
