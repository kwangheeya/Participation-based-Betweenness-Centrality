import time
import sys
from common import get_function_name, debug_print, read_hypergraph, write_bc, adjSimpleEdge
from collections import deque


def bfs_simple(n, s, e):
    # data structure
    dist = {}
    pred = {}
    sigma = {}

    Q = deque()
    S = deque()

    # initialization
    dist[s] = 0;
    sigma[s] = 1;
    pred[s] = set();
    Q.append(s)
    while len(Q) > 0:
        u = Q.popleft()
        S.append(u)
        if u in e:
            for v in e[u]:
                if not (v in dist):
                    dist[v] = dist[u] + 1
                    Q.append(v)
                if dist[v] == dist[u] + 1:
                    if not (v in pred):
                        pred[v] = set()
                    pred[v].add(u)
                    if not (v in sigma):
                        sigma[v] = 0
                    sigma[v] = sigma[v] + sigma[u]
    return pred, sigma, S


def bc_simple(n, e, p_op=False, printFrequency=1000):
    debug_print('START [' + str(get_function_name()) + '] ', print_op=p_op)
    startTime = time.time()
    newE = adjSimpleEdge(n, e)

    BC = [0.0] * n

    for x in range(n):
        value = bc_simple_run(x, n, newE)
        for key in value:
            BC[key] += value[key]/n/(n-1)
        if (x > 0 and x % printFrequency == 0):
            debug_print('>> ' + str(x) + '/' + str(n) + ' ' + str(time.time() - startTime), print_op=p_op)

    endTime = time.time()
    debug_print('> End. Time: ', endTime - startTime, print_op=p_op)
    return BC


def bc_simple_run(node, n, e):
    BC = {}
    pred, sigma, S = bfs_simple(n, node, e)
    delta = {}

    while len(S) > 0:
        u = S.pop()
        if u not in delta:
            delta[u] = 0
        for v in pred[u]:
            if v not in delta:
                delta[v] = 0
            delta[v] = delta[v] + sigma[v] / sigma[u] * (1 + delta[u])
        if u is not node:
            if delta[u] > 0.0:
                BC[u] = delta[u]
    return BC


def bc_simple_main(dataFile='../data/test/test.txt', op='v1', lam=0.5, p_op=True, mode='a'):
    n, e = read_hypergraph(dataFile, mode=mode)
    bc = bc_simple(n, e, p_op=True)
    debug_print(bc, print_op=False)
    temp = (dataFile.split('/')[-1]).split('.')
    # write_bc(bc, "../bc/aBBC_" + temp[0] + '_' + op + '_' + str(lam) + '.' + temp[1])
    write_bc(bc, "../bc/OBC_" + temp[0] + '.' + temp[1])


if __name__ == '__main__':
    # usage: python OBC.py [filenmae]
    # ex) python aBBC.py test.txt v1 0.5
    # ex) python aBBC.py test.txt v2
    inCommand = sys.argv[1:]
    if (len(inCommand) >= 1):
        bc_simple_main(dataFile=inCommand[0])
    else:
        bc_simple_main()