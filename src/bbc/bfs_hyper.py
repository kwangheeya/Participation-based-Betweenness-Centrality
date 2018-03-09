from collections import deque


def bfs_hyper(x, e, fstar):
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
    Q.append(x)
    while len(Q) > 0:
        v = Q.popleft()
        if v in fstar:
            for fIndex in fstar[v]:
                w = e[fIndex][1]
                if not (w in dist):
                    dist[w] = dist[v] + 1
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


def bfs_hyper2(x, e, fstar):
    dist = {}
    bstar = {}
    preh = {}
    sigmaN = {}
    sigmaE = {}
    Q = deque()
    Se = set()  # visited hyperedges
    Sv = deque()  # visited nodes

    # initialization
    dist[x] = 0
    sigmaN[x] = 1
    bstar[x] = set()
    Q.append(x)
    while len(Q) > 0:
        v = Q.popleft()
        Sv.append(v)
        if v in fstar:
            for fIndex in fstar[v]:
                w = e[fIndex][1]
                if not (w in dist):
                    dist[w] = dist[v] + 1
                    Q.append(w)
                if dist[w] == dist[v] + 1:
                    Se.add(fIndex)
                    sigmaE[fIndex] = sigmaE.get(fIndex, 0) + sigmaN[v]
                    sigmaN[w] = sigmaN.get(w, 0) + sigmaN[v]
                    if not (w in bstar):
                        bstar[w] = set()
                    bstar[w].add(fIndex)
                    if not (fIndex in preh):
                        preh[fIndex] = set()
                    preh[fIndex] |= bstar[v]
    return preh, Se, sigmaN, sigmaE
