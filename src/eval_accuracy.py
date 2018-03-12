from compute_bc import get_bcpath
from bbc import BBC
import statistics


def precisionAtN(true_list, sample_list, n):
    true_index_list = sorted(range(len(true_list)), key=true_list.__getitem__, reverse=True)[:n]
    sample_index_list = sorted(range(len(sample_list)), key=sample_list.__getitem__, reverse=True)[:n]
    return len(set(true_index_list).intersection(sample_index_list)) / n


def compute_accuracy():
    filepaths = ["data/HLMN/hlmn.txt", "data/PID/pid.txt", "data/DBLP/dblp5.txt", "data/DBLP/dblp10.txt",
                 "data/CODA/kegg.txt"]
    for filepath in filepaths[:3]:
        bbc = BBC()
        bbc.load_graph(filepath)
        n = bbc.graph.n

        bbc.load_bc(get_bcpath(filepath, 'bbc'))
        bbc_list = bbc.bc[:]
        bbc.load_bc(get_bcpath(filepath, 'lbbc'))
        lbbc_list = bbc.bc[:]
        bbc.load_bc(get_bcpath(filepath, 'slbbc', eps=0.01))
        slbbc_list = bbc.bc[:]

        diff_bbc_lbbc = [abs(x - y) for x, y in zip(bbc_list, lbbc_list)]
        diff_bbc_slbbc = [abs(x - y) for x, y in zip(bbc_list, slbbc_list)]
        diff_lbbc_slbbc = [abs(x-y) for x,y in zip(lbbc_list, slbbc_list)]
        print('BBC-LBBC: max ', max(diff_bbc_lbbc), ", avg ", sum(diff_bbc_lbbc) / n,
              statistics.stdev(diff_bbc_lbbc))
        print('BBC-SLBBC: max ', max(diff_bbc_slbbc), ", avg ", sum(diff_bbc_slbbc) / n,
              statistics.stdev(diff_bbc_slbbc))
        print('LBBC-SLBBC: max ', max(diff_lbbc_slbbc), ", avg ", sum(diff_lbbc_slbbc) / n,
              statistics.stdev(diff_lbbc_slbbc))

        for rn in [1, 2, 3, 5, 10, 50]:
            print('--')
            print('lBBC\t Precision @', rn, ':', precisionAtN(bbc_list, lbbc_list, rn))
            print('slBBC\t Precision @', rn, ':', precisionAtN(bbc_list, slbbc_list, rn))

if __name__ == '__main__':
    compute_accuracy()
