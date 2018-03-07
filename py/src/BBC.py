from collections import deque
from common import get_function_name, debug_print, read_hypergraph, write_bc
from bfs_hyper import bfs_hyper, bfs_hyper2
import time
import sys
import copy


# b-hypergraph
def bc_hyper_single(n, e, op='exact', p_op=False, printFrequency=1000, lam=0):
    # n:|V|, e:hyperedges
    # op: option='exact' or 'apprx'
    if op != 'exact' and op != 'apprx':
        debug_print('Option is wrong')
        return
    else:
        debug_print('START [' + str(get_function_name()) + '] ' + str(op) + ' ' + str(printFrequency) + ' ' + str(p_op),
                    lam, print_op=p_op)

    startTime = time.time()

    fstar = {};

    for i, edge in enumerate(e):
        for v in edge[0]:
            if not (v in fstar):
                fstar[v] = set()
            fstar[v].add(i)

    bc = [0.0 for x in range(n)]

    for v in range(n):
        preh, dist, Se, sigmaN, sigmaE = bfs_hyper(v, e, fstar)

        setSe = set(Se)
        for u in range(n):
            A = setSe & fstar[u] if u in fstar else set()
            A = A - fstar[v] if v in fstar else A
            if len(A) == 0:
                continue;
            S = Se.copy()
            delta = {}
            while len(S) > 0:
                currentIndex = S.pop()
                currentTarget = e[currentIndex][1]
                if not (currentIndex in delta):
                    delta[currentIndex] = 0.0
                delta[currentIndex] += float(sigmaE[currentIndex]) / sigmaN[currentTarget]
                if currentIndex not in A:
                    for preIndex in preh[currentIndex]:
                        if not (preIndex in delta):
                            delta[preIndex] = 0.0
                        delta[preIndex] += float(delta[currentIndex] * sigmaE[preIndex]) / sigmaE[currentIndex]
                else:
                    bc[u] += delta[currentIndex]/n/(n-1)

        if (v > 0 and v % printFrequency == 0):
            debug_print('>> ' + str(v) + '/' + str(n) + ' ' + str(time.time() - startTime), print_op=p_op)
    endTime = time.time()
    debug_print('> End. Time: ' + str(endTime - startTime), print_op=p_op)
    return bc



def bc_hyper_single_main(dataFile='../data/test/test.txt', op='exact', lam=0.5, p_op=True, mode='a'):
    n, e = read_hypergraph(dataFile, mode=mode)
    bc = bc_hyper_single(n, e, op=op, lam=lam, p_op=p_op)
    debug_print(bc, print_op=False)
    temp = (dataFile.split('/')[-1]).split('.')
    #write_bc(bc, "../bc/aBBC_" + temp[0] + '_' + op + '_' + str(lam) + '.' + temp[1])
    write_bc(bc, "../bc/BBC_" + temp[0] +  '.' + temp[1])


if __name__ == '__main__':
    # usage: python aBBC.py [filenmae] [option] [lambda]
    # ex) python aBBC.py test.txt v1 0.5
    # ex) python aBBC.py test.txt v2
    inCommand = sys.argv[1:]
    if (len(inCommand) >= 1):
        bc_hyper_single_main(dataFile=inCommand[0])
    else:
        bc_hyper_single_main()
