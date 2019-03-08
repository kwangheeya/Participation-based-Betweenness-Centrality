from .utils import fopen, print_status
from .graph import Bhypergraph, Graph
from collections import deque
from scipy import optimize
import time
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
        print('\t>Load bc, {0} elements'.format(len(self.bc)))

    def store_bc(self, filepath):
        with fopen(filepath, mode='w') as file:
            for val in self.bc:
                file.write(str(val) + ',')
        print('\t>Store bc, {0} elements'.format(len(self.bc)))

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
        start_time = time.time()
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
                print_status("\t[" + str(time.time()-start_time) + "s]--" + str(x) + '/' + str(n), print_op=print_op)
        if print_op:
            print_status('\t>Total ' + '%.3f' % (time.time() - start_time) + 's', print_op=print_op)

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


class PBC(BC):
    def __init__(self, graph=None, bcoption='pbc', is_get_bc_hyperedge=False):
        super().__init__(graph)
        self.bcoption = bcoption
        self.is_get_bc_hyperedge = is_get_bc_hyperedge

    def store_bcE(self, filepath):
        assert self.bcE, 'No bc scores of hyperedges'
        with fopen(filepath, mode='w') as file:
            for val in self.bcE:
                file.write(str(val) + ',')
        print('\t>Store bc, {0} elements'.format(len(self.bcE)))

    def load_bcE(self, filepath):
        assert self.graph, 'No graph'
        with fopen(filepath) as file:
            line = file.readline()
            line = line[:-1].split(',')
            idx = 0
            lene = len(self.graph.e)
            for text in line:
                if idx >= lene:
                    print('\t-The BC file is wrong')
                    idx = 0
                    self.bcE = None
                    break
                self.bcE[idx] = float(text)
                idx += 1
        print('\t>Load bc, {0} elements'.format(len(self.bcE)))

    def load_graph(self, filepath):
        self.graph = Bhypergraph(filepath)
        self.bc = [0.0] * self.graph.n
        if self.is_get_bc_hyperedge:
            self.bcE = [0.0] * len(self.graph.e)
        else:
            self.bcE = None

    def compute(self, print_op=True, possible_nodes=[], eps=0.01, eta=0.1):
        assert self.graph, 'No graph'
        if self.bcoption.lower() == 'slpbc':
            print_status('+ Start computing {0} with eps={1} eta={2}'.format(self.bcoption.upper(), str(eps), str(eta)),
                         print_op=print_op)
        elif self.bcoption.lower() == 'naiveslpbc':
            print_status('+ Start computing {0} with eps={1} eta={2}'.format(self.bcoption.upper(), str(eps), str(eta)),
                         print_op=print_op)
        elif self.bcoption.lower() == 'bbc':
            print_status('+ Start computing {0} with lambda={1}'.format(self.bcoption.upper(), str(eps)),
                         print_op=print_op)
        else:
            print_status('+ Start computing {0}'.format(self.bcoption.upper()), print_op=print_op)
        n = self.graph.n
        maxval = n//100 if n > 100 else n
        start_time = None
        if print_op:
            start_time = time.time()

        if self.bcoption.lower() == 'pbc':
            self._compute_pbc(start_time, print_op=print_op)
        elif self.bcoption.lower() == 'lpbc':
            self._compute_lpbc(start_time, print_op=print_op)
        elif self.bcoption.lower() == 'hpbc':
            self._compute_hpbc(start_time, print_op=print_op)
        elif self.bcoption.lower() == 'slpbc':
            self._compute_slpbc(start_time, print_op=print_op, possible_nodes=possible_nodes, eps=eps, eta=eta)
        elif self.bcoption.lower() == 'naiveslpbc':
            self._compute_naive_slpbc(start_time, print_op=print_op, possible_nodes=possible_nodes, eps=eps, eta=eta)

        elif self.bcoption.lower() == 'bbc':
            self._compute_bbc(start_time, print_op=print_op, lam=eps)

        if print_op:
            print_status('\t>Total ' + '%.3f' % (time.time() - start_time) + 's', print_op=print_op)

    def _compute_bbc(self, start_time, print_op=True, lam=0.5):
        n = self.graph.n
        maxval = int(n*0.1)
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar
        for v in range(n):
            preh, dist, Se, sigmaN, sigmaE = self._bfs(v)
            delta = {}
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
                            self.bc[u] += lam * delta[currentIndex] / n / (n - 1)
            if print_op and (v % (maxval+1) == 0):
                print_status("\t[" + str(time.time() - start_time) + "s]--" + str(v) + '/' + str(n), print_op=print_op)
    '''
    def _compute_naive_slpbc(self, start_time, print_op=True, possible_nodes=[], eps=0.01, eta=0.1):
        n = self.graph.n
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar
        bstar = self.graph.bstar

        total_num_samples = n * (n - 1)
        prev_num_samples = 0
        s1 = np.int_((1 + 8 * eps + np.sqrt(1 + 16 * eps)) * np.log(6 / eta) / 4 / np.square(eps))
        next_num_samples = s1
        pd = int(np.log2(n))+1
        if next_num_samples > total_num_samples:
            print_status('Too small for sampling. Turn lpbc', print_op=print_op)
            self._compute_lpbc(start_time, print_op=print_op)
            return
        bfs_time = 0
        w_time = 0
        bcE_l2_sqaure = {}
        i = 1
        while True:
            print(i, next_num_samples, '\t>' + '%.3f' % (time.time() - start_time) + 's')
            for s in range(next_num_samples-prev_num_samples):
                [x] = np.random.randint(n, size=1)
                [y] = np.random.randint(n, size=1)
                while x == y:
                    [y] = np.random.randint(n, size=1)
                temp_time = time.time()
                preh, dist, Se, sigmaN, sigmaE = self._bfs_btw_two_nodes(x, y)
                temp_time = time.time() - temp_time
                bfs_time += temp_time
                if y in dist:
                    sigma = {}
                    Rx = set(Se)
                    newSe = deque()
                    e_record = set()
                    for bstar_edge in Rx & bstar[y]:
                        newSe.append(bstar_edge)
                        e_record.add(bstar_edge)

                    while len(newSe) > 0:
                        currentIndex = newSe.pop()
                        currentTarget = e[currentIndex][1]
                        if not (currentIndex in sigma):
                            sigma[currentIndex] = 0.0
                        if currentTarget == y:
                            sigma[currentIndex] = sigmaE[currentIndex]
                        #sigma[currentIndex] = sigmaE[currentIndex] * (1 + sigma[currentIndex])
                        if currentIndex not in bcE_l2_sqaure:
                            bcE_l2_sqaure[currentIndex] = 0.0
                        bcE_l2_sqaure[currentIndex] += (sigma[currentIndex]/sigmaN[y]) * (sigma[currentIndex]/sigmaN[y])
                        if self.is_get_bc_hyperedge:
                            self.bcE[currentIndex] += sigma[currentIndex]/sigmaN[y]
                        for preIndex in preh[currentIndex]:
                            if preIndex not in e_record:
                                newSe.append(preIndex)
                                e_record.add(preIndex)
                            if not (preIndex in sigma):
                                sigma[preIndex] = 0.0
                            sigma[preIndex] += sigma[currentIndex] * sigmaE[preIndex] / sigmaE[currentIndex]

                        if x not in e[currentIndex][0]:
                            d = dist[e[currentIndex][1]]
                            delta = sigma[currentIndex]/sigmaN[y]
                            for u in e[currentIndex][0]:
                                if u in dist:
                                    if d == dist[u] + 1:
                                        self.bc[u] += delta
                                else:
                                    lam_f = 0
                                    for _e in fstar[u] - fstar[x]:
                                        if _e in Rx:
                                            lam_f += 1
                                    lam_f = 1 / lam_f
                                    self.bc[u] += lam_f * delta
            if print_op:
                print_status('\t>cummulatvie bfs time ' + '%.3f' % (bfs_time) + 's', print_op=print_op)
            cal_sign = True
            for j in range(i):
                w_j = 1/(eps*eps)*(pd+j*np.log(2)+np.log(3/eta))
                if (w_j >= s1*((1.5)**(j-1))) and (w_j <= s1*((1.5)**(j))):
                    cal_sign = False
                    break
            if cal_sign:
                temp_time = time.time()
                w_optimal = self._find_w_optimal(bcE_l2_sqaure, next_num_samples)
                temp_time = time.time() - temp_time
                w_time += temp_time
                gamm = np.log(3 / eta) + i * np.log(2)
                xi = 2 * w_optimal + (gamm + np.sqrt(
                    gamm * (gamm + 4 * next_num_samples * w_optimal))) / next_num_samples + np.sqrt(
                    gamm / 2 / next_num_samples)

                if xi < eps:
                    break
                else:
                    prev_num_samples = next_num_samples
                    next_num_samples = PBC.next_sample_size(next_num_samples)
                    i += 1
                    if next_num_samples > n*(n-1):
                        next_num_samples = n*(n-1)
            else:
                break
        if print_op:
            print_status('\t>cummulatvie w time ' + '%.3f' % (w_time) + 's', print_op=print_op)
        print_status('\t>Total sample ', next_num_samples, print_op=print_op)
        for v in range(n):
            self.bc[v] = self.bc[v] / next_num_samples
        if self.is_get_bc_hyperedge:
            for key in range(len(e)):
                self.bcE[key] = self.bcE[key] / next_num_samples
    '''
    def _compute_naive_slpbc(self, start_time, print_op=True, possible_nodes=[], eps=0.01, eta=0.1):
        n = self.graph.n
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar
        bstar = self.graph.bstar

        total_num_samples = n * (n - 1)
        prev_num_samples = 0
        s1 = np.int_((1 + 8 * eps + np.sqrt(1 + 16 * eps)) * np.log(6 / eta) / 4 / np.square(eps))
        next_num_samples = s1
        pd = int(np.log2(n))+1
        if next_num_samples > total_num_samples:
            print_status('Too small for sampling. Turn lpbc', print_op=print_op)
            self._compute_lpbc(start_time, print_op=print_op)
            return
        bfs_time = 0
        w_time = 0
        bcE_l2_sqaure = {}
        i = 1
        while True:
            print(i, next_num_samples, '\t>' + '%.3f' % (time.time() - start_time) + 's')
            for s in range(next_num_samples-prev_num_samples):
                [x] = np.random.randint(n, size=1)
                [y] = np.random.randint(n, size=1)
                while x == y:
                    [y] = np.random.randint(n, size=1)
                temp_time = time.time()
                preh, dist, Se, sigmaN, sigmaE = self._bfs_btw_two_nodes(x, y)
                temp_time = time.time() - temp_time
                bfs_time += temp_time
                if y in dist:
                    delta = {}
                    Rx = set(Se)
                    e_record = set()
                    while len(Se) > 0:
                        currentIndex = Se.pop()
                        currentTarget = e[currentIndex][1]
                        if dist[currentTarget] >= dist[y] and currentTarget is not y:
                            continue
                        if not (currentIndex in delta):
                            delta[currentIndex] = 0.0
                        if currentTarget == y:
                            delta[currentIndex] = sigmaE[currentIndex]/sigmaN[y]
                            if delta[currentIndex] > 1:
                                print('wrong delta')

                        if not (currentIndex in bcE_l2_sqaure):
                            bcE_l2_sqaure[currentIndex] = 0.0
                        bcE_l2_sqaure[currentIndex] += delta[currentIndex] * delta[currentIndex]
                        if self.is_get_bc_hyperedge:
                            self.bcE[currentIndex] += delta[currentIndex]
                        for preIndex in preh[currentIndex]:
                            if not (preIndex in delta):
                                delta[preIndex] = 0.0
                            delta[preIndex] += delta[currentIndex] * sigmaE[preIndex] / sigmaE[currentIndex]

                        if x not in e[currentIndex][0]:
                            d = dist[e[currentIndex][1]]
                            for u in e[currentIndex][0]:
                                if u in dist:
                                    if d == dist[u] + 1:
                                        self.bc[u] += delta[currentIndex]
                                else:
                                    lam_f = 0
                                    for _e in fstar[u] - fstar[x]:
                                        if _e in Rx:
                                            lam_f += 1
                                    lam_f = 1 / lam_f
                                    self.bc[u] += lam_f * delta[currentIndex]
            if print_op:
                print_status('\t>cummulatvie bfs time ' + '%.3f' % (bfs_time) + 's', print_op=print_op)
            cal_sign = True
            for j in range(i):
                w_j = 1/(eps*eps)*(pd+j*np.log(2)+np.log(3/eta))
                if (w_j >= s1*((1.5)**(j-1))) and (w_j <= s1*((1.5)**(j))):
                    cal_sign = False
                    break
            if cal_sign:
                temp_time = time.time()
                w_optimal = self._find_w_optimal(bcE_l2_sqaure, next_num_samples)
                temp_time = time.time() - temp_time
                w_time += temp_time
                gamm = np.log(3 / eta) + i * np.log(2)
                xi = 2 * w_optimal + (gamm + np.sqrt(
                    gamm * (gamm + 4 * next_num_samples * w_optimal))) / next_num_samples + np.sqrt(
                    gamm / 2 / next_num_samples)

                if xi < eps:
                    break
                else:
                    prev_num_samples = next_num_samples
                    next_num_samples = PBC.next_sample_size(next_num_samples)
                    i += 1
                    if next_num_samples > n*(n-1):
                        next_num_samples = n*(n-1)
            else:
                break
        if print_op:
            print_status('\t>cummulatvie w time ' + '%.3f' % (w_time) + 's', print_op=print_op)
        print_status('\t>Total sample ', next_num_samples, print_op=print_op)
        for v in range(n):
            self.bc[v] = self.bc[v] / next_num_samples
        if self.is_get_bc_hyperedge:
            for key in range(len(e)):
                self.bcE[key] = self.bcE[key] / next_num_samples

    def _find_w_optimal(self, bcE_l2_sqaure, next_num_samples):
        w = 0
        t = len(self.graph.e) - len(bcE_l2_sqaure)
        def f(x):
            formula = t
            for v2_key in bcE_l2_sqaure:
                v2 = bcE_l2_sqaure[v2_key]
                formula += np.exp((x*x)*v2/2/(next_num_samples*next_num_samples))
            formula = 1/x * np.log(formula)
            return formula
        result = optimize.minimize(f, np.array([10]), method='L-BFGS-B', options={'disp': None, 'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-08, 'maxfun': 100, 'maxiter': 10000, 'iprint': -1, 'maxls': 20})
        if result.x[0] <= 0:
            print('wrong1')
        w = f(result.x[0])
        return w

    def _compute_slpbc(self, start_time, print_op=True, possible_nodes=[], eps=0.01, eta=0.1):
        n = self.graph.n
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar

        total_num_samples = n * (n - 1)
        prev_num_samples = 0
        s1 = np.int_((1 + 8 * eps + np.sqrt(1 + 16 * eps)) * np.log(6 / eta) / 4 / np.square(eps))
        next_num_samples = s1
        if next_num_samples > total_num_samples:
            print_status('Too small for sampling. Turn lpbc', print_op=print_op)
            self._compute_lpbc(start_time, print_op=print_op)
            return
        bfs_time = 0
        w_time = 0

        VS = {}
        num_possible_nodes = len(possible_nodes) if len(possible_nodes) > 0 else n
        selected_pairs = {}
        i = 1
        while True:
            l = 0
            while l + prev_num_samples < next_num_samples and l + prev_num_samples < total_num_samples:
                [idx_x] = np.random.randint(num_possible_nodes, size=1)
                x = PBC.remapping_node_id(possible_nodes, idx_x)
                if x not in selected_pairs:
                    selected_pairs[x] = {x}
                while len(selected_pairs[x]) == num_possible_nodes:
                    [idx_x] = np.random.randint(num_possible_nodes, size=1)
                    x = PBC.remapping_node_id(possible_nodes, idx_x)

                [idx_y] = np.random.randint(num_possible_nodes, size=1)
                y = PBC.remapping_node_id(possible_nodes, idx_y)
                while y in selected_pairs[x]:
                    [idx_y] = np.random.randint(num_possible_nodes, size=1)
                    y = PBC.remapping_node_id(possible_nodes, idx_y)
                temp_time = time.time()
                preh, dist, Se, sigmaN, sigmaE = self._bfs_btw_two_nodes(x, y)
                temp_time = time.time() - temp_time
                bfs_time += temp_time
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
            print_status('\t>cummulatvie bfs time ' + '%.3f' % (bfs_time) + 's', print_op=print_op)
            if l + prev_num_samples > next_num_samples:
                next_num_samples = l + prev_num_samples
            temp_time = time.time()
            w_optimal = np.sqrt(max(VS.values()) * 2 * np.log(len(e))) / next_num_samples
            temp_time = time.time() - temp_time
            w_time += temp_time
            gamm = np.log(3 / eta) + i * np.log(2)
            delta_s = 2 * w_optimal + (gamm + np.sqrt(
                gamm * (gamm + 4 * next_num_samples * w_optimal))) / next_num_samples + np.sqrt(
                gamm / 2 / next_num_samples)

            if delta_s < eps:
                break
            elif next_num_samples >= total_num_samples:
                break
            else:
                i += 1
                prev_num_samples = next_num_samples
                next_num_samples = PBC.next_sample_size(next_num_samples)


        print_status('\t>cummulatvie w time ' + '%.3f' % (w_time) + 's', print_op=print_op)
        print_status('\t>Total sample ', next_num_samples, print_op=print_op)
        for i in range(n):
            self.bc[i] = self.bc[i] / next_num_samples

        if self.is_get_bc_hyperedge:
            for key in VS:
                self.bcE[key] = VS[key] / next_num_samples

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

    def _compute_pbc(self, start_time, print_op=True):
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
                    if (currentIndex not in A):
                        for preIndex in preh[currentIndex]:
                            if not (preIndex in delta):
                                delta[preIndex] = 0.0
                            delta[preIndex] += float(delta[currentIndex] * sigmaE[preIndex]) / sigmaE[currentIndex]
                    else:
                        self.bc[u] += delta[currentIndex] / n / (n - 1)
            if print_op and (v % (maxval+1) == 0):
                print_status("\t[" + str(time.time() - start_time) + "s]--" + str(v) + '/' + str(n), print_op=print_op)

    def _compute_lpbc(self, start_time, print_op=True):
        n = self.graph.n
        maxval = n // 100 if n > 100 else n
        self.bc = [0.0] * n
        e = self.graph.e
        fstar = self.graph.fstar
        for v in range(n):
            preh, dist, Se, sigmaN, sigmaE = self._bfs(v)
            delta = {}
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
                                if (_e in preh) and len(preh[_e]) > 0:
                                    l += 1
                            l = 1 / l
                            self.bc[u] += l * delta[currentIndex] / n / (n - 1)
            if print_op and (v % (maxval+1) == 0):
                print_status("\t[" + str(time.time() - start_time) + "s]--" + str(v) + '/' + str(n), print_op=print_op)

    def _compute_hpbc(self, start_time, print_op=True):
        n = self.graph.n
        maxval = int(n*0.1)
        e = self.graph.e
        self.bc = [0.0] * len(e)
        fstar = self.graph.fstar
        for v in range(n):
            preh, dist, Se, sigmaN, sigmaE = self._bfs(v)
            delta = {}
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
                self.bc[currentIndex] += delta[currentIndex] / n / (n - 1)
            if print_op and (v % (maxval+1) == 0):
                print_status("\t[" + str(time.time() - start_time) + "s]--" + str(v) + '/' + str(n), print_op=print_op)

    def compute_from_seed(self, v, print_op=True):
        assert self.graph, 'No graph'
        assert v < self.graph.n, 'Wrong seed'
        print_status('+ Start computing LPBC from seed node {0}'.format(str(v)), print_op=print_op)
        n = self.graph.n
        maxval = n//100 if n > 100 else n

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
