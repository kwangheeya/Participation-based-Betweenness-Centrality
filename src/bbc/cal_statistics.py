from collections import deque
from common import read_hypergraph, read_bc
import math
import statistics

def covariance(X, Y):
    if len(X) == 0 or len(Y) == 0:
        return None
    mean_X = math.fsum(X) / len(X)
    mean_Y = math.fsum(Y) / len(Y)
    cov = 0
    for i in range(len(X)):
        cov += (X[i] - mean_X) * (Y[i] - mean_Y)
    cov = cov / len(X)
    return cov


def spearman(rankX, rankY):
    cov = covariance(rankX, rankY)
    devX = math.sqrt(covariance(rankX, rankX))
    devY = math.sqrt(covariance(rankY, rankY))
    r = cov / (devX * devY)
    return r


def precisionAtN(true_list, sample_list, n):
    true_index_list = sorted(range(len(true_list)), key=true_list.__getitem__, reverse=True)[:n]
    sample_index_list = sorted(range(len(sample_list)), key=sample_list.__getitem__, reverse=True)[:n]
    #print(true_index_list, sample_index_list)
    return len(set(true_index_list).intersection(sample_index_list)) / (n)


def getRank(l):
    n = len(l)
    indices = list(range(n))
    indices.sort(key=lambda x: -l[x])
    rank = [0] * n
    for i, x in enumerate(indices):
        rank[x] = i
    return rank


def statistics_main():
    '''
    dataFile = "pueraria.txt"
    n, e = read_hypergraph("data/" + dataFile)
    lams = [0.5, 1.0]
    temp = dataFile.split('.')
    datah = temp[0] + '_' + str(1.0) + '.' + temp[1]
    bc_exact = read_bc("bc/bch_exact_" + datah, n)
    bc_exact_rank = getRank(bc_exact)
    for lam in lams:
        print(lam)
        datah = temp[0] + '_' + str(lam) + '.' + temp[1]
        bc = read_bc("bc/bch_apprx_" + datah, n)

        bc_rank = getRank(bc)
        print('Spearman:', spearman(bc_exact_rank, bc_rank))
        for i in [1, 2, 3, 5, 10]:
            print('Precision @', i, ':', precisionAtN(bc_exact, bc, i))
    '''
    data_files = ["/HLMN/hlmn.txt","/PID/pid.txt","/DBLP/dblp5.txt","/DBLP/dblp10.txt","/CODA/kegg.txt"]
    for data_file in data_files:
        print()
        n, e = read_hypergraph("../data" + data_file)
        startable = set()
        targetable = set()
        for i, edge in enumerate(e):
            for v in edge[0]:
                if not (v in startable):
                    startable.add(v)
            if not (edge[1] in targetable):
                targetable.add(edge[1])
        total_num_pairs = len(startable)*len(targetable)
        temp = (data_file.split('/')[-1]).split(".")
        BBC = read_bc("../bc/BBC_"+temp[0]+"."+temp[1], n)
        #aBBC = read_bc("../bc/aBBC_" + temp[0] + "_v1."+temp[1], n)
        aBBC2 = read_bc("../bc/aBBC_" + temp[0] + "_v2." + temp[1], n)
        if not aBBC2:
            aBBC2 = [0 for x in range(n)]
        sBBC = read_bc("../bc/sBBC_" + temp[0] + "_01."+temp[1], n)

        #diff_BBC_aBBC = [0 for x in range(n)]
        diff_BBC_aBBC2 = [0 for x in range(n)]
        diff_BBC_sBBC = [0 for x in range(n)]
        #diff_aBBC_sBBC = [0 for x in range(n)]
        diff_aBBC2_sBBC = [0 for x in range(n)]
        for i in range(n):
            if BBC is not None:
                #BBC[i] = BBC[i]
                #diff_BBC_aBBC[i] = abs(BBC[i] - aBBC[i])
                diff_BBC_sBBC[i] = abs(BBC[i] - sBBC[i])
                diff_BBC_aBBC2[i] = abs(BBC[i] - aBBC2[i])
            #diff_aBBC_sBBC[i] = abs(sBBC[i] - aBBC[i])
            diff_aBBC2_sBBC[i] = abs(sBBC[i] - aBBC2[i])

        #print('BBC-aBBC: max ', max(diff_BBC_aBBC),", avg ",sum(diff_BBC_aBBC)/n, statistics.stdev(diff_BBC_aBBC))
        print('BBC-aBBC2: max ', max(diff_BBC_aBBC2), ", avg ", sum(diff_BBC_aBBC2) / n, statistics.stdev(diff_BBC_aBBC2))
        print('BBC-sBBC: max ', max(diff_BBC_sBBC), ", avg ", sum(diff_BBC_sBBC) / n, statistics.stdev(diff_BBC_sBBC))
        #print('aBBC-sBBC: max ', max(diff_aBBC_sBBC), ", avg ", sum(diff_aBBC_sBBC) / n)
        print('aBBC2-sBBC: max ', max(diff_aBBC2_sBBC), ", avg ", sum(diff_aBBC2_sBBC) / n, statistics.stdev(diff_aBBC2_sBBC))

        if BBC is not None:
            for rn in [1, 2, 3, 5, 10, 50, 100]:
                print('--')
                #print('aBBC v1\t Precision @', rn, ':', precisionAtN(BBC, aBBC, rn))
                print('aBBC\t Precision @', rn, ':', precisionAtN(BBC, aBBC2, rn))
                print('sBBC\t Precision @', rn, ':', precisionAtN(BBC, sBBC, rn))

if __name__ == '__main__':
    statistics_main()
