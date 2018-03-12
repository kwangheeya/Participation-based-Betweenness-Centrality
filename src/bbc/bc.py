from .utils import fopen, print_status
from .graph import Bhypergraph, Graph
from collections import deque
import progressbar
import numpy as np


class BC:
    def __init__(self, graph=None):
        self.graph = graph
        if graph:
            self.bc = [0.0]*self.graph.n
        else:
            self.bc = None

    def load_bc(self, filepath):
        assert self.graph, 'No graph'
        with fopen(filepath) as file:
            line = file.readline()
            line = line[:-1].split(',')
            idx = 0
            for text in line:
                if idx >= self.graph.n:
                    print('\t-The BC file is wrong')
                    idx = 0
                    self.bc = None
                    break
                self.bc[idx] = float(text)
                idx += 1
        print('\t>Load bc, {0} nodes'.format(len(self.bc)))

    def store_bc(self, filepath):
        with fopen(filepath, mode='w') as file:
            for val in self.bc:
                file.write(str(val) + ',')
        print('\t>Store bc, {0} nodes'.format(len(self.bc)))

    def sorted_index(self, reverse=True):
        assert self.bc, 'No bc list'
        bc = self.bc
        return sorted(range(len(bc)), key=bc.__getitem__, reverse=reverse)

    def load_graph(self, filepath):
        raise NotImplementedError

    def store_graph(self, filepath):
        raise NotImplementedError


class OBC(BC):
    def load_graph(self, filepath):
        self.graph = Graph(filepath)
        self.bc = [0.0] * self.graph.n

    def compute(self, print_op=True, possible_nodes=None, eps=None, eta=None):
        assert self.graph, 'No graph'
        print_status('+ Start computing OBC', print_op=print_op)
        n = self.graph.n
        maxval = n//100 if n > 100 else n
        bar = progressbar.ProgressBar(maxval=maxval+1, widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                                                progressbar.Percentage(), ' ', progressbar.Timer()])
        if print_op:
            bar.start()
        self.bc = [0.0] * n
        for x in range(n):
            pred, sigma, S = self._bfs(x)
            delta = {}
            while len(S) > 0:
                u = S.pop()
                if u not in delta:
                    delta[u] = 0
                for v in pred[u]:
                    if v not in delta:
                        delta[v] = 0
                    delta[v] = delta[v] + sigma[v] / sigma[u] * (1 + delta[u])
                if u is not x:
                    if delta[u] > 0.0:
                        self.bc[u] += delta[u]/n/(n-1)
            if print_op and (x % (maxval+1) == 0):
                bar.update(x//100) if n > 100 else bar.update(x)

        if print_op:
            bar.finish()

    def _bfs(self, s):
        e = self.graph.e
        # data structure
        dist = {}
        pred = {}
        sigma = {}

        Q = deque()
        S = deque()

        # initialization
        dist[s] = 0
        sigma[s] = 1
        pred[s] = set()
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


class BBC(BC):
    def __init__(self, graph=None, bcoption='bbc'):
        super().__init__(graph)
        self.bcoption = bcoption

    def load_graph(self, filepath):
        self.graph = Bhypergraph(filepath)
        self.bc = [0.0] * self.graph.n

    def compute(self, print_op=True, possible_nodes=[], eps=0.01, eta=0.1):
        assert self.graph, 'No graph'
        if self.bcoption.lower() == 'slbbc':
            print_status('+ Start computing {0} with eps={1} eta={2}'.format(self.bcoption.upper(), str(eps), str(eta)),
                         print_op=print_op)
        else:
            print_status('+ Start computing {0}'.format(self.bcoption.upper()), print_op=print_op)
        n = self.graph.n
        maxval = n//100 if n > 100 else n
        bar = progressbar.ProgressBar(max_value=maxval+1, widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                            progressbar.Percentage(), ' ', progressbar.Timer()], redirect_stdout=True)
        if print_op:
            bar.start()

        if self.bcoption.lower() == 'bbc':
            self._compute_bbc(bar, print_op=print_op)
        elif self.bcoption.lower() == 'lbbc':
            self._compute_lbbc(bar, print_op=print_op)
        elif self.bcoption.lower() == 'slbbc':
            self._compute_slbbc(bar, print_op=print_op, possible_nodes=possible_nodes, eps=eps, delta=eta)

        if print_op:
            bar.finish()

    def _compute_slbbc(self, bar, print_op=True, possible_nodes=[], eps=0.01, delta=0.1):
        n = self.graph.n
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar

        total_num_samples = n * (n - 1)
        if len(possible_nodes) > 0:
            total_num_samples = len(possible_nodes) * (len(possible_nodes) - 1)
        next_num_samples = np.int_((1 + 4 * eps + np.sqrt(1 + 8 * eps)) * np.log(6 / delta) / 4 / np.square(eps))

        if len(possible_nodes) >= len(e) or next_num_samples > total_num_samples:
            self._compute_lbbc(bar, print_op=print_op)
        else:
            VS = {}
            prev_num_samples = 0
            num_possible_nodes = len(possible_nodes) if len(possible_nodes) > 0 else n
            selected_pairs = {}
            i = 1
            while True:
                l = 0
                while l + prev_num_samples < next_num_samples and l + prev_num_samples < total_num_samples:
                    [idx_x] = np.random.randint(num_possible_nodes, size=1)
                    x = BBC.remapping_node_id(possible_nodes, idx_x)
                    if x not in selected_pairs:
                        selected_pairs[x] = {x}
                    while len(selected_pairs[x]) == num_possible_nodes:
                        [idx_x] = np.random.randint(num_possible_nodes, size=1)
                        x = BBC.remapping_node_id(possible_nodes, idx_x)

                    [idx_y] = np.random.randint(num_possible_nodes, size=1)
                    y = BBC.remapping_node_id(possible_nodes, idx_y)
                    while y in selected_pairs[x]:
                        [idx_y] = np.random.randint(num_possible_nodes, size=1)
                        y = BBC.remapping_node_id(possible_nodes, idx_y)

                    preh, dist, Se, sigmaN, sigmaE = self._bfs_btw_two_nodes(x, y)

                    if y in dist:
                        delt = {}
                        Rx = set(Se)
                        while len(Se) > 0:
                            currentIndex = Se.pop()
                            currentTarget = e[currentIndex][1]
                            if not (currentIndex in delt):
                                delt[currentIndex] = 0.0
                            if not (e[currentIndex][1] in selected_pairs[x]):
                                delt[currentIndex] += sigmaE[currentIndex] / sigmaN[currentTarget]
                            if currentIndex not in VS:
                                VS[currentIndex] = 0.0
                            VS[currentIndex] += delt[currentIndex]
                            for preIndex in preh[currentIndex]:
                                if not (preIndex in delt):
                                    delt[preIndex] = 0.0
                                delt[preIndex] += delt[currentIndex] * sigmaE[preIndex] / sigmaE[currentIndex]
                            if x not in e[currentIndex][0]:
                                d = dist[e[currentIndex][1]]
                                for u in e[currentIndex][0]:
                                    if u in dist:
                                        if d == dist[u] + 1:
                                            self.bc[u] += delt[currentIndex]
                                    else:
                                        lam_f = 0
                                        for _e in fstar[u] - fstar[x]:
                                            if _e in Rx:
                                                lam_f += 1
                                        lam_f = 1 / lam_f
                                        self.bc[u] += lam_f * delt[currentIndex]
                    additional_sample = 0
                    if y not in dist:
                        selected_pairs[x].add(y)
                        additional_sample += 1
                    else:
                        for v in dist:
                            if v not in selected_pairs[x]:
                                selected_pairs[x].add(v)
                                additional_sample += 1
                    l = l + additional_sample

                if l + prev_num_samples > next_num_samples:
                    next_num_samples = l + prev_num_samples
                w_optimal = np.sqrt(max(VS.values()) * 2 * np.log(len(e))) / next_num_samples
                gamm = np.log(3 / delta) + i * np.log(2)
                delta_s = 2 * w_optimal + (gamm + np.sqrt(
                    gamm * (gamm + 4 * next_num_samples * w_optimal))) / 2 / next_num_samples + np.sqrt(
                    gamm / 2 / next_num_samples)

                if delta_s < eps:
                    break
                elif next_num_samples >= total_num_samples:
                    break
                else:
                    i += 1
                    prev_num_samples = next_num_samples
                    next_num_samples = BBC.next_sample_size(next_num_samples)

                    if print_op:
                        bar.max_value = next_num_samples // 100
                        bar.update(prev_num_samples // 100)

            print_status('\t>Total sample ', next_num_samples, print_op=print_op)
            for i in range(n):
                self.bc[i] = self.bc[i] / next_num_samples

    @staticmethod
    def remapping_node_id(possible_nodes, x):
        if len(possible_nodes) > x:
            x = possible_nodes[x]
        return x

    @staticmethod
    def next_sample_size(S, c=1.2):
        return int(S * c)

    def _bfs_btw_two_nodes(self, x, y):
        fstar = self.graph.fstar
        e = self.graph.e
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

    def _compute_bbc(self, bar, print_op=True):
        n = self.graph.n
        maxval = n // 100 if n > 100 else n
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar
        for v in range(n):
            preh, dist, Se, sigmaN, sigmaE = self._bfs(v)

            setSe = set(Se)
            for u in range(n):
                A = setSe & fstar[u] if u in fstar else set()
                A = A - fstar[v] if v in fstar else A
                if len(A) == 0:
                    continue
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
                        self.bc[u] += delta[currentIndex] / n / (n - 1)
            if print_op and (v % (maxval+1) == 0):
                bar.update(v//100) if n > 100 else bar.update(v)

    def _compute_lbbc(self, bar, print_op=True):
        n = self.graph.n
        maxval = n // 100 if n > 100 else n
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar
        for v in range(n):
            preh, dist, Se, sigmaN, sigmaE = self._bfs(v)
            delta = {}
            Rx = set(Se)
            while len(Se) > 0:
                currentIndex = Se.pop()
                currentTarget = e[currentIndex][1]
                if not (currentIndex in delta):
                    delta[currentIndex] = 0.0
                delta[currentIndex] += sigmaE[currentIndex] / sigmaN[currentTarget]
                for preIndex in preh[currentIndex]:
                    if not (preIndex in delta):
                        delta[preIndex] = 0.0
                    delta[preIndex] += delta[currentIndex] * sigmaE[preIndex] / sigmaE[currentIndex]
                if v not in e[currentIndex][0]:
                    d = dist[e[currentIndex][1]]
                    for u in e[currentIndex][0]:
                        if u in dist:
                            if d == dist[u] + 1:
                                self.bc[u] += delta[currentIndex] / n / (n - 1)
                        else:
                            l = 0
                            for _e in fstar[u] - fstar[v]:
                                if _e in Rx:
                                    l += 1
                            l = 1 / l
                            self.bc[u] += l * delta[currentIndex] / n / (n - 1)
            if print_op and (v % (maxval+1) == 0):
                bar.update(v//100) if n > 100 else bar.update(v)

    def compute_from_seed(self, v, print_op=True):
        assert self.graph, 'No graph'
        assert v < self.graph.n, 'Wrong seed'
        print_status('+ Start computing LBBC from seed node {0}'.format(str(v)), print_op=print_op)
        n = self.graph.n
        maxval = n//100 if n > 100 else n
        bar = progressbar.ProgressBar(max_value=maxval+1, widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                            progressbar.Percentage(), ' ', progressbar.Timer()], redirect_stdout=True)
        if print_op:
            bar.start()

        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar

        preh, dist, Se, sigmaN, sigmaE = self._bfs(v)
        delta = {}
        Rx = set(Se)
        while len(Se) > 0:
            currentIndex = Se.pop()
            currentTarget = e[currentIndex][1]
            if not (currentIndex in delta):
                delta[currentIndex] = 0.0
            delta[currentIndex] += sigmaE[currentIndex] / sigmaN[currentTarget]
            for preIndex in preh[currentIndex]:
                if not (preIndex in delta):
                    delta[preIndex] = 0.0
                delta[preIndex] += delta[currentIndex] * sigmaE[preIndex] / sigmaE[currentIndex]
            if v not in e[currentIndex][0]:
                d = dist[e[currentIndex][1]]
                for u in e[currentIndex][0]:
                    if u in dist:
                        if d == dist[u] + 1:
                            self.bc[u] += delta[currentIndex] / n / (n - 1)
                    else:
                        l = 0
                        for _e in fstar[u] - fstar[v]:
                            if _e in Rx:
                                l += 1
                        l = 1 / l
                        self.bc[u] += l * delta[currentIndex] / n / (n - 1)

        if print_op:
            bar.finish()

    def _bfs(self, x):
        fstar = self.graph.fstar
        e = self.graph.e
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
