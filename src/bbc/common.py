from collections import deque, Iterable
import inspect
import logging
import traceback
import os.path

def debug_print(*obj, print_op=True):
    if print_op & isinstance(obj, Iterable):
        for i in obj:
            print(i, end=' ')
        print()


def get_function_name():
    return traceback.extract_stack(None, 2)[0][2]


def read_hypergraph(_filename, mode='b'):
    HYPERGRAPH_file = str(_filename)
    try:
        HYPERGRAPH_fp = open(HYPERGRAPH_file, 'r');
    except IOError:
        print('Error: can\'t find ', HYPERGRAPH_file, ' or read data')
        return None
    else:
        print('Reading data', '[' + HYPERGRAPH_file + ']')

    e = []

    line = HYPERGRAPH_fp.readline()
    n = int(line[:-1])

    for line in HYPERGRAPH_fp:
        temp = line[:-1].split(';')
        sourceT = temp[0].split(',')
        targetT = temp[1].split(',')
        source = set()
        for text in sourceT:
            source.add(int(text))
        if mode == 'b':
            for text in targetT:
                target = int(text)
                if len(source) > 0 and target not in source:
                    e.append((source, target))
        elif mode == 'a':
            for text in targetT:
                target = int(text)
                if len(source) > 0 and not (target in source and len(source) ==1):
                    e.append((source, target))
        else:
            target = set()
            for text in targetT:
                target.add(int(text))
            if len(source) > 0 and len(target) > 0 and source != target:
                if (source, target) not in e:
                    e.append((source, target))
                    # else:
                    #	print((source,target))

    print('> Done. #nodes:', n, ' #hyperedges :', len(e))
    HYPERGRAPH_fp.close()
    return n, e


def rw_hypergraph(_filename, mode='h'):
    HYPERGRAPH_file = './' + str(_filename)
    hyperedge_file = './data/dblp5.txt'
    try:
        HYPERGRAPH_fp = open(HYPERGRAPH_file, 'r');
    except IOError:
        print('Error: can\'t find ', HYPERGRAPH_file, ' or read data')
        return None
    else:
        print('Reading data', '[' + HYPERGRAPH_file + ']')

    e = ()

    line = HYPERGRAPH_fp.readline()
    n = int(line[:-1])

    for line in HYPERGRAPH_fp:
        temp = line[:-1].split(';')
        sourceT = temp[0].split(',')
        targetT = temp[1].split(',')
        source = set()
        for text in sourceT:
            source.add(int(text))
        if mode == 'b':
            for text in targetT:
                target = int(text)
                if len(source) > 0 and target not in source:
                    if (source, target) not in e:
                        e = e + ((source, target),)
        else:
            target = set()
            for text in targetT:
                target.add(int(text))
            if len(source) > 0 and len(target) > 0 and source != target:
                if (source, target) not in e:
                    e = e + ((source, target),)
                    # else:
                    #	print((source,target))

    print('> Done. #nodes:', n, ' #hyperedges :', len(e))
    HYPERGRAPH_fp.close()

    try:
        hyperedge_fp = open(hyperedge_file, 'w');
    except IOError:
        print('Error: can\'t find ', hyperedge_file, ' or write data')
        return None
    else:
        print(hyperedge_file)
    hyperedge_fp.write(str(n) + '\n')
    for edge in e:
        src = ''
        for line in edge[1]:
            src += str(line) + ','
        for line in edge[0]:
            hyperedge_fp.write(str(src[:-1]) + ';' + str(line) + '\n')
    hyperedge_fp.close()


def adjSimpleEdge(n, e, mode='b'):
    newE = {}
    if mode == 'b':
        for edge in e:
            for u in edge[0]:
                if not (u in newE):
                    newE[u] = set()
                v = edge[1]
                if u != v:
                    newE[u].add(v)
    else:
        for edge in e:
            for u in edge[0]:
                if not (u in newE):
                    newE[u] = set()
                for v in edge[1]:
                    if u != v:
                        newE[u].add(v)
    return newE


def adjUSimpleEdge(n, e):
    newE = [0] * n
    for i in range(n):
        newE[i] = set()
    for edge in e:
        for u in edge[0]:
            v = edge[1]
            if u != v:
                newE[u].add(v)
                newE[v].add(u)
    return newE


def findCC(s, n, e):
    # data structure
    visitedN = [False] * n
    Q = deque()
    S = deque()

    # initialization
    visitedN[s] = True;
    Q.append(s)
    while len(Q) > 0:
        # print Q
        u = Q.popleft()
        S.append(u)
        for v in e[u]:
            if not visitedN[v]:
                visitedN[v] = True
                Q.append(v)
    return S


def lcc_size(n, e):
    largest_component = 0
    S = set()
    newE = adjUSimpleEdge(n, e)
    for s in range(n):
        if s not in S:
            T = findCC(s, n, newE)
            if len(T) > largest_component:
                largest_component = len(T)
            # print largest_component
            S |= set(T)
    return largest_component


def read_bc(_filename, n):
    BC_file = str(_filename)
    try:
        BC_fp = open(BC_file, 'r');
    except IOError:
        print('Error: can\'t find ', BC_file, ' or read data')
        return None

    bc = []
    line = BC_fp.readline()
    line = line[:-1].split(',')
    for text in line:
        bc.append(float(text))

    BC_fp.close()

    if len(bc) != n:
        print('Error : Read bc')
        return -1

    print('Read bc ', '[' + BC_file + '],', len(bc), 'nodes')

    return bc


def write_bc(bc, _filename):
    BC_file = str(_filename)
    BC_fp = open(BC_file, 'w');

    for val in bc:
        BC_fp.write(str(val) + ',')

    BC_fp.close()

    print('Written bc', '[' + BC_file + '],', len(bc), 'nodes')


def mean_source(e):
    s = 0
    maxS = 1
    minS = None
    for edge in e:
        cardi = len(edge[0])
        s += cardi
        minS = cardi if minS is None else min(minS, cardi)
        maxS = max(maxS, cardi)
    return minS, maxS, s / len(e)


def read_entity2id(_file, mode='text'):
    _file = './' + str(_file)
    _fp = open(_file, 'r');
    _fp.readline()
    entity2id = dict()
    for line in _fp:
        temp = line[:-1].split('\t')
        key = temp[0]
        value = int(temp[1])
        if mode != 'text':
            key = int(key)
        entity2id[key] = value
    return entity2id


def find_entity_name(global_entity2id, local_entity2id, value):
    intermediate_value = None
    final_key = ''
    for key in local_entity2id:
        if local_entity2id[key] == value:
            intermediate_value = key
            break;
    for key in global_entity2id:
        if global_entity2id[key] == intermediate_value:
            final_key = key
            break
    return final_key

def dblp_hanlding():
    data_file = "/DBLP/dblp5p.txt"
    n, e = read_hypergraph("../data" + data_file, mode='a')

    fstar = {}
    bstar = {}

    for i, edge in enumerate(e):
        for v in edge[0]:
            if not (v in fstar):
                fstar[v] = set()
            fstar[v].add(i)
        if edge[1] not in bstar:
            bstar[edge[1]] = set()
        bstar[edge[1]].add(i)

    for v in fstar:
        if len(fstar[v]) > 300 and v not in bstar:
            print(v, len(fstar[v]))

    data_file = "/DBLP/dblp5.txt"
    n, e = read_hypergraph("../data" + data_file, mode='a')
    cnt_bbc = 0
    cnt_obc = 0
    desired_node = 2490
    for edge in e:
        if desired_node in edge[0]:
            cnt_bbc += 1
        elif desired_node == edge[1]:
            cnt_obc += 1
    print(cnt_bbc, cnt_obc)

def hlmn_hanlding():
    data_file = "/HLMN/hlmn.txt"
    n, e = read_hypergraph("../data" + data_file, mode='a')
    temp = (data_file.split('/')[-1]).split(".")

    k= 245

    bbc = read_bc("../bc/BBC_"+temp[0]+"."+temp[1], n)
    sorted_index = sorted(range(len(bbc)), key=lambda i: -bbc[i])
    top20bbc = set(sorted_index[:k])
    obc = read_bc("../bc/OBC_"+temp[0]+"."+temp[1], n)
    sorted_index = sorted(range(len(obc)), key=lambda i: -obc[i])
    top20obc = set(sorted_index[:k])
    cnt_bbc = 0
    cnt_obc = 0
    for edge in e:
        if len(edge[0] & top20bbc) > 0 or edge[1] in top20bbc:
            cnt_bbc += 1
        if len(edge[0] & top20obc) > 0 or edge[1] in top20obc:
            cnt_obc += 1
    print(cnt_bbc, cnt_obc, len(top20bbc & top20obc))

def common_main():
    dataFiles = ["/HLMN/hlmn.txt","/PID/pid.txt","/DBLP/dblp5.txt","/DBLP/dblp10.txt","/CODA/kegg.txt","/bitcoin/bitcoin600k.txt"]
    # rw_hypergraph("data/"+dataFile)

    for dataFile in dataFiles:
        n, e = read_hypergraph("../data" + dataFile, mode='a')
        # rw_hypergraph("data/"+dataFile)
        #print('lcc_size(n, e)', lcc_size(n, e))
        print('reduction', sum([len(x) for x in adjSimpleEdge(n, e).values()]))

        s = 0
        max_src = 0
        for edge in e:
            if len(edge[0]) > 1:
                s += 1
            if len(edge[0]) > max_src:
                max_src = len(edge[0])
        print('NON-SINGLETON hyperedges', s)
        print('MAX_SRC ', max_src)
    #print(mean_source(e))
    '''
    collectionOfSourceSets = ();

    for i, edge in enumerate(e):
        if not (edge[0] in collectionOfSourceSets):
            collectionOfSourceSets += (edge[0],)
    print('#Source set ', len(collectionOfSourceSets))
    

    dataFile = "pueraria.txt"
    global_entity2id = read_entity2id("data/entity2id.txt")
    local_entity2id = read_entity2id("data/entity2id_" + dataFile, mode='int')
    n, e = read_hypergraph("data/" + dataFile)

	indegree= [0]*n
	outdegree = [0]*n
	for edge in e:
		for u in edge[0]:
			outdegree[u] += 1
		indegree[edge[1]] += 1
	BC_file = './'+str("degree.txt")
	BC_fp = open(BC_file, 'w');

	for val in indegree:
		BC_fp.write(str(val)+',')
	BC_fp.write('\n')
	for val in outdegree:
		BC_fp.write(str(val)+',')	
	BC_fp.close()

    lams = [0.5, 1.0]
    temp = dataFile.split('.')
    rank = []
    for lam in lams:
        datah = temp[0] + '_' + str(lam) + '.' + temp[1]
        bc = read_bc("bc/bch_apprx_" + datah, n)
        top = sorted(range(len(bc)), key=lambda i: bc[i])
        gene = set()

        while (len(gene) < 16):
            value = top.pop()
            name = find_entity_name(global_entity2id, local_entity2id, value)
            if name.find('GE') > -1:
                gene.add(name[:10])
                print(name)
                # gene.add(name[12:-1])#

        rank.append(gene)
    rank = rank[0] & rank[1]
    print(len(rank), rank)

    
    n, e = read_hypergraph("data/" + dataFile)
    lams = [0.5, 1.0]
    rank = []
    temp = dataFile.split('.')
    for lam in lams:
        datah = temp[0] + '_' + str(lam) + '.' + temp[1]
        bc = read_bc("bc/bch_apprx_" + datah, n)
        top = sorted(range(len(bc)), key=lambda i: -bc[i])
        rank.append(top)
    for i in range(int(n * 0.1)):
        if rank[0][i] != rank[1][i]:
            print(i)
    '''


if __name__ == '__main__':
    #dblp_hanlding()
    hlmn_hanlding()
    #common_main()
