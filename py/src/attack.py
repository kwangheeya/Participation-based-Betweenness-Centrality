from aBBC import bc_hyper_single
from OBC import bc_simple
from sBBC import sBBC_compute
from common import get_function_name, debug_print, read_hypergraph, read_bc, lcc_size
import sys


def write_result(u, val, lcc, method, dataFile):
    CC_file = ''
    temp = dataFile.split('.')
    if method == 'r':
        CC_file = '../result/v_' + temp[0] + '_' + method + '.' + temp[1]
    elif method == 's':
        temp2 = temp[0].split('_')
        CC_file = '../result/v_' + temp2[0] + '_' +  method + '_' + temp2[1] + '.' + temp[1]
    else:
        CC_file = '../result/v_' + temp[0] + '.' + temp[1]

    CC_fp = open(CC_file, 'a');
    CC_fp.write(str(u) + "\t" + str(val) + "\t" + str(lcc) + "\n")
    CC_fp.close()


def delete_node(n, e, bc):
    maxIndex = bc.index(max(bc))
    newE = ()
    for edge in e:
        if not (maxIndex in edge[0]) and maxIndex != edge[1]:
            newE += (edge,)
    n_inodes = 0

    non_inodes = set()
    for edge in newE:
        non_inodes |= edge[0]
        non_inodes.add(edge[1])
    n_inodes = n - len(non_inodes)

    return maxIndex, bc[maxIndex], n_inodes, newE, list(non_inodes)


def attack_vulnerability(n, e, dataFile, method='r', lam=0, printFrequency=50):
    debug_print('START [' + str(get_function_name()) + '] ', method, lam)
    if method == 'v1' or method =='v2':
        temp = (dataFile.split('/')[-1]).split('.')
        dataFile = temp[0] + '_' + method + '.' + temp[1]
    elif method == 'r':
        dataFile = dataFile.split('/')[-1]
    elif method == 's':
        temp = (dataFile.split('/')[-1]).split('.')
        dataFile = temp[0] + '_' + str(lam).split('.')[-1] + '.' + temp[1]
    write_result('index', 'bc_value', 'lcc\tinodes', method, dataFile)
    val = 1
    for i in range(n):
        if i % printFrequency == 0:
            debug_print('>', i)
        bc = []

        if method == 'r':
            bc = bc_simple(n, e) if i > 0 else read_bc("../bc/OBC_" + dataFile, n)
        elif method == 'v1' or method == 'v2':
            bc = bc_hyper_single(n, e, op=method) if i > 0 else read_bc("../bc/aBBC_" + dataFile, n)
        elif method == 's':
            '''
            if val + lam < 1:
                eps = lam * (val + lam)
            else:
                eps = lam
            '''
            eps = lam
            bc = sBBC_compute(n, e, possible_nodes=possible_nodes, eps=eps) if i > 0 else read_bc("../bc/sBBC_" + dataFile, n)

        u, val, n_inodes, e, possible_nodes = delete_node(n, e, bc)
        # print (bc_list,e)
        # print e
        lcc = lcc_size(n, e)
        write_result(u, val, str(lcc) + '\t' + str(n_inodes), method, dataFile)
        if lcc == 1 or val == 0:
            break
    debug_print('END [' + str(get_function_name()) + '] ')


def attack_main(dataFile='../data/test/test.txt', method='r', lams=[]):
    for lam in reversed(lams):
        n, e = read_hypergraph(dataFile)
        attack_vulnerability(n, e, dataFile, method=method, lam=lam)


if __name__ == '__main__':
    # usage: python attack.py [filenmae] [option]
    # ex) python bc_hyper_single.py test.txt v1
    inCommand = sys.argv[1:]
    if len(inCommand) >= 2:
        if inCommand[1] == 'r':
            attack_main(dataFile=inCommand[0], method=inCommand[1], lams=[0])
        elif inCommand[1] == 'v1' or inCommand[1] == 'v2':
            attack_main(dataFile=inCommand[0], method=inCommand[1], lams=[0])
        elif inCommand[1] == 's' and len(inCommand) >= 3:
            print(float(inCommand[2]))
            attack_main(dataFile=inCommand[0], method=inCommand[1], lams=[float(inCommand[2])])
    else:
        attack_main()
