from collections import deque
from common import get_function_name, debug_print, read_hypergraph, write_bc
from aBBC import bc_hyper_single
import time
import sys
import copy
import numpy as np


def bfs_for_sBBC(x, y, e, fstar):
    dist = {}
    bstar = {}
    preh = {}
    sigmaN = {}
    sigmaE = {}
    Q = deque()
    S = deque()  # stack for hyperedges
    visited = set()
    # initialization
    dist[x] = 0
    sigmaN[x] = 1
    bstar[x] = set()
    if x != y:
        Q.append(x)
    exit_d = -1
    while len(Q) > 0 and exit_d:
        v = Q.popleft()
        if v in fstar:
            for fIndex in fstar[v]:
                w = e[fIndex][1]
                if not (w in dist):
                    dist[w] = dist[v] + 1
                    if w == y:
                        exit_d = dist[w]
                    if exit_d > 0 and dist[w] > exit_d:
                        exit_d = False
                        break
                    bstar[w] = set()
                    Q.append(w)
                if dist[w] == dist[v] + 1:
                    if not (fIndex in visited):
                        visited.add(fIndex)
                        S.append(fIndex)
                    sigmaE[fIndex] = sigmaE.get(fIndex, 0) + sigmaN[v]
                    sigmaN[w] = sigmaN.get(w, 0) + sigmaN[v]
                    bstar[w].add(fIndex)
                    if not (fIndex in preh):
                        preh[fIndex] = set()
                    preh[fIndex] |= bstar[v]
    return preh, dist, S, sigmaN, sigmaE

def sBBC_compute(n, e, eps=1, delta=0.1, op='v2', lam=0.5, p_op=False, printFrequency=1000):
    debug_print('START [' + str(get_function_name()) + '] ', print_op=p_op)
    debug_print('\tOption : ' + str(eps) + ' ' + str(delta) + ' ' + str(printFrequency) + ' ' + str(p_op), print_op=p_op)

    startTime = time.time()
    fstar = {}

    startable = set()
    targetable = set()
    for i, edge in enumerate(e):
        for v in edge[0]:
            if not (v in fstar):
                fstar[v] = set()
                startable.add(v)
            fstar[v].add(i)
        if not (edge[1] in targetable):
            targetable.add(edge[1])
    startable = list(startable)
    targetable = list(targetable)

    bc = [0.0 for x in range(n)]

    total_num_samples = len(startable)*len(targetable)
    next_num_samples = np.int_((1 + 4 * eps + np.sqrt(1 + 8 * eps)) * np.log(6 / delta) / 4 / np.square(eps))
    print('\t next_num_samples ',next_num_samples,' / total ',total_num_samples)
    if next_num_samples >= total_num_samples * 0.5:
        print('Too small for sampling. Turn to aBBC')
        bc = bc_hyper_single(n, e, p_op=p_op, lam=lam, op=op)
    else:
        VS = [0.0 for x in range(len(e))]
        prev_num_samples = 0
        selected_pairs = {}

        i = 1
        while(True):
            debug_print('\tsample', next_num_samples, print_op=p_op)
            l = 0
            while l+prev_num_samples < next_num_samples and l+prev_num_samples < total_num_samples * 0.5:
                [idx] = np.random.randint(len(startable), size = 1)
                [idy] = np.random.randint(len(targetable), size=1)
                x = startable[idx]
                y = targetable[idy]

                resample_count = 0
                while resample_count < 1000 and x in selected_pairs and y in selected_pairs[x]:
                    [idx] = np.random.randint(len(startable), size=1)
                    [idy] = np.random.randint(len(targetable), size=1)
                    x = startable[idx]
                    y = targetable[idy]
                    resample_count+= 1

                if x not in selected_pairs:
                    selected_pairs[x] = set()

                if resample_count == 1000:
                    print('wow')

                preh, dist, Se, sigmaN, sigmaE = bfs_for_sBBC(x, y, e, fstar)

                if x!=y and y in dist:
                    delt = {}
                    Rx = set(Se)
                    while len(Se) > 0:
                        currentIndex = Se.pop()
                        currentTarget = e[currentIndex][1]
                        if not (currentIndex in delt):
                            delt[currentIndex] = 0.0
                        if not (e[currentIndex][1] in selected_pairs[x]):
                            delt[currentIndex] += float(sigmaE[currentIndex]) / sigmaN[currentTarget]
                        VS[currentIndex] += delt[currentIndex]
                        for preIndex in preh[currentIndex]:
                            if not (preIndex in delt):
                                delt[preIndex] = 0.0
                            delt[preIndex] += float(delt[currentIndex] * sigmaE[preIndex]) / sigmaE[currentIndex]
                        if x not in e[currentIndex][0]:
                            d = dist[e[currentIndex][1]]
                            for u in e[currentIndex][0]:
                                if u in dist:
                                    if d == dist[u] + 1:
                                        bc[u] += delt[currentIndex]
                                else:
                                    lam_f = lam
                                    if op == 'v2':
                                        _l = 0
                                        for _e in fstar[u] - fstar[x]:
                                            if _e in Rx:
                                                _l += 1
                                        lam_f = 2 / (_l + 1)
                                    bc[u] += lam_f * delt[currentIndex]
                additional_sample = 0
                if (x ==y) or (y not in dist):
                    selected_pairs[x].add(y)
                    additional_sample += 1
                else:
                    for v in dist:
                        if v != x and v not in selected_pairs[x]:
                            selected_pairs[x].add(v)
                            additional_sample += 1
                l = l + additional_sample

            if l+prev_num_samples > next_num_samples:
                next_num_samples = l+prev_num_samples
            w_optimal = np.sqrt(max(VS)*2*np.log(len(e)))/next_num_samples
            gamm = np.log(3 / delta) + i*np.log(2)
            delta_s = 2*w_optimal + (gamm+np.sqrt(gamm*(gamm+4*next_num_samples*w_optimal)))/2/next_num_samples+np.sqrt(gamm/2/next_num_samples)
            debug_print('\toptimal_w',w_optimal,'\tdetla_s', delta_s, print_op=p_op)

            if delta_s < eps:
                break
            elif next_num_samples >= total_num_samples:
                break
            else:
                i += 1
                prev_num_samples = next_num_samples
                next_num_samples = next_sample_size(next_num_samples)
        print('\tTotal sample ', next_num_samples)
        for i in range(n):
            bc[i] = bc[i] / next_num_samples

    endTime = time.time()
    debug_print('> End. Time: ' + str(endTime - startTime), print_op=p_op)
    return bc

def next_sample_size(S, c=1.2):
    return int(S*c)


def sBBC_main(dataFile='../data/test/test.txt', eps=[0.01], delta=0.1, p_op=True, mode='a'):
    n, e = read_hypergraph(dataFile, mode=mode)
    for epsilon in eps:
        bc = sBBC_compute(n, e, epsilon, delta, p_op=p_op)
        debug_print(bc, print_op=False)
        temp = (dataFile.split('/')[-1]).split('.')
        write_bc(bc, "../bc/sBBC_" + temp[0] + '_' + str(epsilon).split('.')[-1] + '.' + temp[1])


if __name__ == '__main__':
    # usage: python sBBC.py [filenmae] [epsilons] [delta]
    # ex) python sBBC.py test.txt 0.2
    # ex) python sBBC.py test.txt 0.2 0.1
    # ex) python sBBC.py test.txt 0.2,0.21 0.1
    inCommand = sys.argv[1:]
    if len(inCommand) >= 2:
        eps_str = inCommand[1].split(',')
        eps = []
        for ep_str in eps_str:
            ep = float(ep_str)
            eps.append(ep)
        if len(inCommand) >= 3:
            sBBC_main(dataFile=inCommand[0], eps=eps, delta=float(inCommand[2]))
        else:
            sBBC_main(dataFile=inCommand[0], eps=eps)
    else:
        sBBC_main()
