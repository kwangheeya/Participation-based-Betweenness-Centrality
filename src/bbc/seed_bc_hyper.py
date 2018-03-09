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

    seed_set = read_seedset('data/coda2_herbal_compound_idmap.txt')

    for v in seed_set:
        if(v>=n):
            continue;
        preh, dist, Se, sigmaN, sigmaE = bfs_hyper(v, e, fstar)
        '''
		setSe = set(Se)
		for u in range(n):
			A = setSe & fstar[u] if u in fstar else set()
			A = A - fstar[v] if v in fstar else A
			if len(A) == 0:
				continue;
			S = Se.copy()
			delta = {}
			while len(S)>0:
				currentIndex = S.pop()
				currentTarget = e[currentIndex][1]
				if not (currentIndex in delta):
					delta[currentIndex] = 0.0
				delta[currentIndex] += float(sigmaE[currentIndex])/sigmaN[currentTarget]
				if currentIndex not in A:
					for preIndex in preh[currentIndex]:
						if not (preIndex in delta):
							delta[preIndex] = 0.0				
						delta[preIndex] += float(delta[currentIndex]*sigmaE[preIndex])/sigmaE[currentIndex]
			for currentIndex in A:
				if v not in e[currentIndex][0]:
					bc[u] += delta[currentIndex]
		'''
        delta = {}
        while len(Se) > 0:
            currentIndex = Se.pop()
            currentTarget = e[currentIndex][1]
            if not (currentIndex in delta):
                delta[currentIndex] = 0.0
            delta[currentIndex] += float(sigmaE[currentIndex]) / sigmaN[currentTarget]
            for preIndex in preh[currentIndex]:
                if not (preIndex in delta):
                    delta[preIndex] = 0.0
                delta[preIndex] += float(delta[currentIndex] * sigmaE[preIndex]) / sigmaE[currentIndex]
            if v not in e[currentIndex][0]:
                d = dist[e[currentIndex][1]]
                for u in e[currentIndex][0]:
                    if (u in dist) and d == dist[u] + 1:
                        bc[u] += delta[currentIndex]
                    if not (u in dist):
                        bc[u] += lam * delta[currentIndex]

        # print(v,delta)
        # print(v,sigmaE)
        # print(v,sigmaN)
        # debug_print('>>2 ', time.time()-startTime,print_op=p_op)
        if (v > 0 and v % printFrequency == 0):
            debug_print('>> ' + str(v) + '/' + str(n) + ' ' + str(time.time() - startTime), print_op=p_op)
    endTime = time.time()
    debug_print('> End. Time: ' + str(endTime - startTime), print_op=p_op)
    return bc

def read_seedset(_filename):
    HYPERGRAPH_file = './' + str(_filename)
    try:
        HYPERGRAPH_fp = open(HYPERGRAPH_file, 'r');
    except IOError:
        print('Error: can\'t find ', HYPERGRAPH_file, ' or read data')
        return None
    else:
        print('Reading data', '[' + HYPERGRAPH_file + ']')

    seed_set = set()

    for line in HYPERGRAPH_fp:
        temp = line[:-1].split('\t')
        seed_id = temp[0]
        seed_set.add(int(seed_id))

    print('> Done. #seed nodes:', len(seed_set))
    HYPERGRAPH_fp.close()
    return seed_set


def bc_hyper_single_main(dataFile='test.txt', op='apprx', lam=0, p_op=True, mode='b'):
    n, e = read_hypergraph("data/" + dataFile, mode=mode)
    bc = bc_hyper_single(n, e, op=op, lam=lam, p_op=p_op)
    debug_print(bc, print_op=False)
    temp = dataFile.split('.')
    write_bc(bc, "bc/bch_" + op + '_' + temp[0] + '_' + str(lam) + '.' + temp[1])


if __name__ == '__main__':
    # usage: python bc_hyper_single.py [filenmae] [lambda] [option]
    # ex) python bc_hyper.py test.txt 0.5 a
    inCommand = sys.argv[1:]
    if (len(inCommand) >= 3):
        temp = inCommand[1].split(',')
        for s in temp:
            lam = float(s)  # if inCommand[1].isnumeric() else 0.5
            bc_hyper_single_main(dataFile=inCommand[0], lam=lam, mode=inCommand[2])
    elif (len(inCommand) >= 2):
        temp = inCommand[1].split(',')
        for s in temp:
            lam = float(s)  # if inCommand[1].isnumeric() else 0.5
            bc_hyper_single_main(dataFile=inCommand[0], lam=lam)
    elif (len(inCommand) >= 1):
        bc_hyper_single_main(dataFile=inCommand[0])
    else:
        bc_hyper_single_main()
