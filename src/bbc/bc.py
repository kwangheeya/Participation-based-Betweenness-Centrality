from .utils import fopen, print_status
from .graph import Bhypergraph, Graph
from collections import deque
import progressbar
from time import sleep, time


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

    def sort_bc(self, order=''):
        assert self.bc, 'No bc list'
        pass

    def load_graph(self, filepath):
        pass

    def store_graph(self, filepath):
        pass


class OBC(BC):
    def load_graph(self, filepath):
        self.graph = Graph(filepath)
        self.bc = [0.0] * self.graph.n

    def remove_edge(self, edgelist):
        pass

    def compute(self, print_op=True):
        assert self.graph, 'No graph'
        print_status('+ Start computing OBC', print_op=print_op)
        n = self.graph.n
        maxval = n//100 if n > 100 else n
        bar = progressbar.ProgressBar(maxval=maxval+1, widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                                                progressbar.Percentage(), ' ', progressbar.ETA()])
        if print_op:
            bar.start()
        self.bc = [0.0] * n
        for x in range(n):
            pred, sigma, S = self.bfs(x)
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
            if print_op & x % (maxval+1) == 0:
                bar.update(x//100) if n > 100 else bar.update(x)

        if print_op:
            bar.finish()

    def bfs(self, s):
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
    pass
