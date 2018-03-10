from collections import deque
from .utils import fopen


class Graph:
    def __init__(self, filepath=None):
        self.n = 0
        self.e = {}
        if filepath:
            self.load(filepath)

    def load(self, filepath):
        with fopen(filepath) as file:
            line = file.readline()
            self.n = int(line[:-1])
            for line in file:
                temp = line[:-1].split(';')
                sourceT = temp[0].split(',')
                targetT = temp[1].split(',')
                source = set()
                for sT in sourceT:
                    source = int(sT)
                    if source not in self.e:
                        self.e[source] = set()
                    for tT in targetT:
                        target = int(tT)
                        if source != target:
                            self.e[source].add(target)
        print('\t>Load #nodes: {0}, #edges : {1}'.format(self.n, sum([len(self.e[x]) for x in self.e])))

    def store(self):
        pass

    def make_undir(self):
        e = {}
        for u in self.e:
            for v in self.e[u]:
                if u not in e:
                    e[u] = set()
                if v not in e:
                    e[v] = set()
                e[u].add(v)
                e[v].add(u)
        return e

    def lcc(self):
        largest_component = 0
        S = set()
        e = self.make_undir()
        for s in range(self.n):
            if s not in S:
                T = self.find_cc(s, e)
                if len(T) > largest_component:
                    largest_component = len(T)
                S |= T
        return largest_component

    def find_cc(self, s, e):
        visited = set()
        Q = deque()
        visited.add(s)
        Q.append(s)
        while len(Q) > 0:
            u = Q.popleft()
            if u in e:
                for v in e[u]:
                    if v not in visited:
                        visited.add(v)
                        Q.append(v)
        return visited


class Bhypergraph(Graph):
    def __init__(self, filepath=None):
        self.n = 0
        self.e = []
        self.fstar = None
        if filepath:
            self.load(filepath)

    def make_undir(self):
        e = {}
        for edge in self.e:
            v = edge[1]
            for u in edge[0]:
                if u not in e:
                    e[u] = set()
                if v not in e:
                    e[v] = set()
                e[u].add(v)
                e[v].add(u)
        return e

    def load(self, filepath):
        with fopen(filepath) as file:
            line = file.readline()
            self.n = int(line[:-1])
            for line in file:
                temp = line[:-1].split(';')
                sourceT = temp[0].split(',')
                targetT = temp[1].split(',')
                source = set()
                for text in sourceT:
                    source.add(int(text))
                for text in targetT:
                    target = int(text)
                    if len(source) > 0 and not (target in source and len(source) == 1):
                        self.e.append((source, target))
        print('\t>Load #nodes: {0}, #hyperedges : {1}'.format(self.n, len(self.e)))

