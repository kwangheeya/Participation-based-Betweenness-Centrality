from common import get_function_name, debug_print, read_hypergraph, read_bc
import sys
import time

def write_result(u, val, lcc, method, dataFile):
    CC_file = ''
    temp = dataFile.split('.')
    if method == 'r':
        CC_file = '../result/r_' + temp[0] + '_' + method + '.' + temp[1]
    elif method == 's':
        temp2 = temp[0].split('_')
        CC_file = '../result/r_' + temp2[0] + '_' +  method + '_' + temp2[1] + '.' + temp[1]
    elif method == 'v1' or method == 'v2':
        CC_file = '../result/r_' + temp[0] + '.' + temp[1]

    CC_fp = open(CC_file, 'a')
    CC_fp.write(str(u) + "\t" + str(val) + "\t" + str(lcc) + "\n")
    CC_fp.close()


def bfs_hyperpath(x,A,e,fstar):
    reachable_nodes = set([x])
    prev_n = len(reachable_nodes)

    remaining_hyperedges = set([edgeindex for edgeindex in range(len(e))])
    if x in fstar:
        for edgeindex in fstar[x]:
            d_set = e[edgeindex][0] - A - reachable_nodes
            if len(d_set) == 0:
                remaining_hyperedges.remove(edgeindex)
                reachable_nodes.add(e[edgeindex][1])
    next_n = len(reachable_nodes)

    u_set = A | reachable_nodes
    while(prev_n < next_n):
        prev_n = len(reachable_nodes)
        beremoval_hyperedges = set()
        for edgeindex in remaining_hyperedges:
            if e[edgeindex][1] in reachable_nodes:
                beremoval_hyperedges.add(edgeindex)
                continue
            if e[edgeindex][0].issubset(u_set) and len(e[edgeindex][0] & reachable_nodes) > 0:
                beremoval_hyperedges.add(edgeindex)
                reachable_nodes.add(e[edgeindex][1])
                u_set.add(e[edgeindex][1])
        next_n = len(reachable_nodes)
        remaining_hyperedges = remaining_hyperedges - beremoval_hyperedges
    return reachable_nodes


def countReachable(n, e, dataFile, method='r', lam=0, printFrequency=50, op=2):
    debug_print('START [' + str(get_function_name()) + '] ', method, lam)
    startTime = time.time()
    bc = []
    if method == 'v1' or method =='v2':
        temp = (dataFile.split('/')[-1]).split('.')
        dataFile = temp[0] + '_' + method + '.' + temp[1]
        bc = read_bc("../bc/aBBC_" + dataFile, n)
    elif method == 'r':
        dataFile = dataFile.split('/')[-1]
        bc = read_bc("../bc/OBC_" + dataFile, n)
    elif method == 's':
        temp = (dataFile.split('/')[-1]).split('.')
        dataFile = temp[0] + '_' + str(lam).split('.')[-1] + '.' + temp[1]
        bc = read_bc("../bc/sBBC_" + dataFile, n)

    write_result('k', 'Satisfiable nodes', '', method, dataFile)


    no_bstar_nodes = set([v for v in range(n)])
    fstar = {}
    for i, edge in enumerate(e):
        for v in edge[0]:
            if not (v in fstar):
                fstar[v] = set()
            fstar[v].add(i)
        if edge[1] in no_bstar_nodes:
            no_bstar_nodes.remove(edge[1])

    if op == 1:
        default_num_reachable_nodes = [0 for x in range(n)]
        for x in range(n):
            T_A = bfs_hyperpath(x, set(), e, fstar)
            default_num_reachable_nodes[x] = len(T_A)

        debug_print('\t compute for each node ' + str(time.time() - startTime))

        sorted_index = sorted(range(len(bc)), key=lambda i: -bc[i])
        for i in range(10):
            end = int(len(sorted_index) * 0.01 * (i + 1))
            A = set(sorted_index[:end])
            num_addi_reachable_nodes = 0
            for x in range(n):
                T_A = bfs_hyperpath(x, A, e, fstar)
                num_addi_reachable_nodes += (len(T_A) - default_num_reachable_nodes[x])/n
            write_result((i+1)*1, num_addi_reachable_nodes, '', method, dataFile)
            debug_print('\t top '+str((i+1)*1)+'% nodes ' + str(time.time() - startTime))
    else:
        print('\t #no_bstar_nodes', len(no_bstar_nodes))
        default_num_reachable_nodes = {}
        i=0
        for x in no_bstar_nodes:
            i += 1
            T_A = bfs_hyperpath(x, set(), e, fstar)
            default_num_reachable_nodes[x] = len(T_A)
            if i % 10000 == 0:
                debug_print('\t\t  '+str(i)+',' + str(time.time() - startTime))
        print(sum([default_num_reachable_nodes[x] for x in default_num_reachable_nodes])/len(no_bstar_nodes))
        debug_print('\t compute for each node ' + str(time.time() - startTime))

        sorted_index = sorted(range(len(bc)), key=lambda i: -bc[i])

        for i in range(10):
            end = int(len(sorted_index) * 0.03 * (i + 1))
            A = set(sorted_index[:end])
            num_addi_reachable_nodes = 0
            for x in no_bstar_nodes:
                T_A = bfs_hyperpath(x, A, e, fstar)
                num_addi_reachable_nodes += (len(T_A) - default_num_reachable_nodes[x]) / len(no_bstar_nodes)
            write_result((i+1)*3, num_addi_reachable_nodes, '', method, dataFile)
            debug_print('\t top ' + str((i+1)*3) + '% nodes ' + str(time.time() - startTime))
    endTime = time.time()
    debug_print('> End. Time: ' + str(endTime - startTime))

def reachability_main(dataFile='../data/test/test.txt', method='r', lams=[]):
    for lam in reversed(lams):
        n, e = read_hypergraph(dataFile)
        countReachable(n, e, dataFile, method=method, lam=lam)


if __name__ == '__main__':
    # usage: python attack.py [filenmae] [option]
    # ex) python reachability.py ../data/test/test.txt v2
    inCommand = sys.argv[1:]
    if len(inCommand) >= 2:
        if inCommand[1] == 'r':
            reachability_main(dataFile=inCommand[0], method=inCommand[1], lams=[0])
        elif inCommand[1] == 'v1' or inCommand[1] == 'v2':
            reachability_main(dataFile=inCommand[0], method=inCommand[1], lams=[0])
        elif inCommand[1] == 's' and len(inCommand) >= 3:
            print(float(inCommand[2]))
            reachability_main(dataFile=inCommand[0], method=inCommand[1], lams=[float(inCommand[2])])
    else:
        reachability_main()
