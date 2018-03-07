from common import get_function_name, debug_print, read_hypergraph, read_bc
import sys
import time

def write_result(u, val, lcc, method, dataFile):
    CC_file = ''
    temp = dataFile.split('.')
    if method == 'r':
        CC_file = '../result/s_' + temp[0] + '_' + method + '.' + temp[1]
    elif method == 's':
        temp2 = temp[0].split('_')
        CC_file = '../result/s_' + temp2[0] + '_' +  method + '_' + temp2[1] + '.' + temp[1]
    else:
        CC_file = '../result/s_' + temp[0] + '.' + temp[1]

    CC_fp = open(CC_file, 'a')
    CC_fp.write(str(u) + "\t" + str(val) + "\t" + str(lcc) + "\n")
    CC_fp.close()


def bfs_hyperpath(A,e,fstar):
    prev_n = len(A)
    A_fstar = set()
    for v in A:
        if v in fstar:
            for edgeindex in fstar[v]:
                A_fstar.add(edgeindex)

    remaining_hyperedges = set([edgeindex for edgeindex in range(len(e))])
    for edgeindex in A_fstar:
        d_set = e[edgeindex][0] - A
        if len(d_set) == 0:
            remaining_hyperedges.remove(edgeindex)
            A.add(e[edgeindex][1])
    next_n = len(A)

    while(prev_n < next_n):
        prev_n = len(A)
        beremoval_hyperedges = set()
        for edgeindex in remaining_hyperedges:
            if e[edgeindex][1] in A:
                beremoval_hyperedges.add(edgeindex)
                continue
            d_set = e[edgeindex][0] - A
            if len(d_set) == 0:
                beremoval_hyperedges.add(edgeindex)
                A.add(e[edgeindex][1])
        next_n = len(A)
        remaining_hyperedges = remaining_hyperedges - beremoval_hyperedges
    return A


def countSatisfiable(n, e, dataFile, method='r', lam=0, printFrequency=50):
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

    #default_nodes = bfs_hyperpath(no_bstar_nodes, e, fstar)
    #default_num_satisfiable_nodes = len(default_nodes)
    default_nodes = set()
    default_num_satisfiable_nodes = 0
    debug_print('\t compute for default ' + str(time.time() - startTime))

    #Compute T_k for each top k% nodes
    sorted_index = sorted(range(len(bc)), key=lambda i: -bc[i])
    for i in range(10):
        end = int(len(sorted_index) * 0.05 * (i + 1))
        A = set(sorted_index[:end])
        default_num_satisfiable_nodes = len(A)
        T_A = bfs_hyperpath(A, e, fstar)
        num_addi_satisfiable_nodes = len(T_A) - default_num_satisfiable_nodes
        write_result(i+1, num_addi_satisfiable_nodes, '', method, dataFile)
        debug_print('\t top '+str(i+1)+'% nodes ' + str(time.time() - startTime))

    endTime = time.time()
    debug_print('> End. Time: ' + str(endTime - startTime))

def countSatisfiable_main(dataFile='../data/test/test.txt', method='r', lams=[]):
    for lam in reversed(lams):
        n, e = read_hypergraph(dataFile)
        countSatisfiable(n, e, dataFile, method=method, lam=lam)


if __name__ == '__main__':
    # usage: python attack.py [filenmae] [option]
    # ex) python countSatisfiable.py ../data/test/test.txt v2
    inCommand = sys.argv[1:]
    if len(inCommand) >= 2:
        if inCommand[1] == 'r':
            countSatisfiable_main(dataFile=inCommand[0], method=inCommand[1], lams=[0])
        elif inCommand[1] == 'v1' or inCommand[1] == 'v2':
            countSatisfiable_main(dataFile=inCommand[0], method=inCommand[1], lams=[0])
        elif inCommand[1] == 's' and len(inCommand) >= 3:
            print(float(inCommand[2]))
            countSatisfiable_main(dataFile=inCommand[0], method=inCommand[1], lams=[float(inCommand[2])])
    else:
        countSatisfiable_main()
