from utils import parseM
from numpy import *
from sys import argv

readCSV = lambda path: [l.split() for l in path.splitlines()]
head = lambda arr: arr[0]
tail = lambda arr: arr[1:]
rptsum = lambda arr: repeat(sum(arr), size(arr))
map2 = lambda fn, mat: map(lambda arr: map(fn, arr), mat)
mapvsum = lambda mat: matrix(map(rptsum, mat))
maphsum = lambda mat: mapvsum(mat).transpose()
idxmin = lambda mat: unravel_index(argmin(mat), shape(mat))
exists = lambda x: x != None
wraparr = lambda x: [x]
# Remove rows and columns of matrix with the listed indices
withoutIndices = lambda m, ids: delete(delete(m, ids, axis=0), ids, axis=1)
# Append a vector as both a row and a column
appendRowCol = lambda m, v: hstack((vstack((m, [v])), map(wraparr, v + [0])))

def neighborJoin(D, forest):
    if len(D) == 2:
        #return Tree(forest[0], forest[1])
        return None
    SH = mapvsum(D)
    SV = SH.transpose()
    I = identity(len(D))
    M = D + (multiply(I, SH + SV) - SH - SV) / (len(D) - 2)
    print M
    i, j = idxmin(M)
    u = [(D[i,k] + D[j,k] - D[i,j]) / 2 for k in range(len(D))]
    #su = (D[i,j] + SH[i,j] - SV[i,j]) / 2
    dui = (D[i,j] + SH[i,j] - SV[i,j]) / 2
    duj = (D[j,i] + SH[j,i] - SV[j,i]) / 2
    #forest = hstack((forest, [Tree(forest[i], forest[j], dui, duj)]))
    D = appendRowCol(D, u)
    D = withoutIndices(D, (i, j))
    #print D

    #forest = delete(forest, (i, j))
    return neighborJoin(D, forest)

rlm = parseM('''
0   295 300 524 1077    1080    978 941 940
295 0   314 487 1071    1088    1010    963 966
300 314 0   472 1085    1088    1025    965 956
524 487 472 0   1101    1099    1021    962 965
1076    1070    1085    1101    0   818 1053    1057    1054
1082    1088    1088    1098    818 0   1070    1085    1080
976 1011    1025    1021    1053    1070    0   963 961
941 963 965 962 1057    1085    963 0   16
940 966 956 965 1054    1080    961 16  0''')

#print readCSV(rlm)
wikim=parseM('''
0   5   9   9   8
5   0   10  10  9
9   10  0   8   7
9   10  8   0   3
8   9   7   3   0''')


if __name__ == "__main__":
       #forest = array(map(Leaf, tail(head(cells))))
       forest = None
       D = wikim
       tree = neighborJoin(D, forest)
       print tree
