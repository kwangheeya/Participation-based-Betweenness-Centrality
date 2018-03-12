from collections import deque
from .utils import fopen


class Graph:
    def __init__(self, filepath=None):
        self.n = 0
        self.e = {}
        self.hyperedge = None
        if filepath:
            self.load(filepath)

    def load(self, filepath, load_op='unkeep'):
        if load_op == 'keep':
            self.hyperedge = []
        with fopen(filepath) as file:
            line = file.readline()
            self.n = int(line[:-1])
            for line in file:
                temp = line[:-1].split(';')
                sourceT = temp[0].split(',')
                targetT = temp[1].split(',')
                for sT in sourceT:
                    source = int(sT)
                    for tT in targetT:
                        target = int(tT)
                        if source != target:
                            if source not in self.e:
                                self.e[source] = set()
                            self.e[source].add(target)
                if load_op == 'keep':
                    source = set()
                    for text in sourceT:
                        source.add(int(text))
                    for text in targetT:
                        target = int(text)
                        if len(source) > 0 and not (target in source and len(source) == 1):
                            self.hyperedge.append((source, target))
        print('\t>Load #nodes: {0}, #edges : {1}'.format(self.n, sum([len(self.e[x]) for x in self.e])))

    def store(self):
        pass

    def remove_nodes(self, nid_list):
        if self.hyperedge:
            nid_list = set(nid_list)
            e = []
            for edge in self.hyperedge:
                if (len(edge[0] & nid_list) > 0) or (edge[1] in nid_list):
                    pass
                else:
                    e.append(edge)
            self.hyperedge = e
            e = {}
            for edge in self.hyperedge:
                v = edge[1]
                for u in edge[0]:
                    if u not in e:
                        e[u] = set()
                    e[u].add(v)
            self.e = e
        else:
            self._remove_nodes(nid_list)

    def _remove_nodes(self, nid_list):
        for nid in nid_list:
            if nid in self.e:
                self.e.pop(nid)
        to_be_deleted = []
        for src in self.e:
            self.e[src] = self.e[src] - nid_list
            if len(self.e[src]) == 0:
                to_be_deleted.append(src)
        for src in to_be_deleted:
            self.e.pop(src)

    def find_nodes_having_edges(self):
        non_inodes = set(self.e.keys())
        for s in self.e.values():
            non_inodes |= s
        return non_inodes

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
        self.fstar = {}
        if filepath:
            self.load(filepath)

    def find_nodes_having_edges(self):
        non_inodes = set()
        for edge in self.e:
            non_inodes |= edge[0]
            non_inodes.add(edge[1])
        return non_inodes

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

    def load(self, filepath, load_op='unkeep'):
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
                        for v in source:
                            if not (v in self.fstar):
                                self.fstar[v] = set()
                            self.fstar[v].add(len(self.e)-1)
        print('\t>Load #nodes: {0}, #hyperedges : {1}'.format(self.n, len(self.e)))

    def remove_nodes(self, nid_list):
        nid_list = set(nid_list)
        e = []
        for edge in self.e:
            if (len(edge[0] & nid_list) > 0) or (edge[1] in nid_list):
                pass
            else:
                e.append(edge)
        self.e = e
        self._reconstruct_fstar()

    def _reconstruct_fstar(self):
        e = self.e
        self.fstar = {}
        for i, edge in enumerate(e):
            for v in edge[0]:
                if not (v in self.fstar):
                    self.fstar[v] = set()
                self.fstar[v].add(i)
