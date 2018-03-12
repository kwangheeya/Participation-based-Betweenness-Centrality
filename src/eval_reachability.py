import sys
from bbc import BBC, fopen
from compute_bc import get_bcpath
import progressbar


def get_reachpath(filepath, bcoption, eps=0.01):
    filename = filepath.split('/')[-1]
    _path = 'result/reachability/r_'
    if bcoption.lower() == 'slbbc':
        s_filename = filename.split('.')
        _path += 'SLBBC_' + s_filename[0] + '_' + str(eps).split('.')[-1]+'.'+s_filename[1]
    else:
        _path += bcoption.upper() + '_' + filename
    return _path


def bfs_hyperpath(x, A, e, fstar):
    reachable_nodes = {x}
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


def eval_reachability(filepath, bcoption='slbbc', epsilons=[0.01], etas=[0.1]):
    bc = BBC()
    bc.load_graph(filepath)
    if bcoption == 'slbbc':
        for eps in epsilons:
            for eta in etas:
                bc_path = get_bcpath(filepath, bcoption, eps=eps)
                bc.load_bc(bc_path)
                bc_sorted_idx = bc.sorted_index()
                reachpath = get_reachpath(filepath, bcoption, eps=eps)
                print('+ Start computing reachability in {0} with eps={1} eta={2}'.format(bcoption.upper(), eps, eta))
                _compute_reachability(bc.graph, bc_sorted_idx, reachpath)
    else:
        bc_path = get_bcpath(filepath, bcoption)
        bc.load_bc(bc_path)
        bc_sorted_idx = bc.sorted_index()
        reachpath = get_reachpath(filepath, bcoption)
        print('+ Start computing reachability in {0}'.format(bcoption.upper()))
        _compute_reachability(bc.graph, bc_sorted_idx, reachpath)


def _compute_reachability(graph, bc_sorted_idx, reachpath, start=0, mode='w'):
    n = graph.n
    e = graph.e
    fstar = graph.fstar
    bstar_nodes = set()
    for edge in e:
        bstar_nodes.add(edge[1])
    no_bstar_nodes = set([v for v in range(n) if v not in bstar_nodes])
    print('\t>#no_bstar_nodes', len(no_bstar_nodes))
    bar = progressbar.ProgressBar(maxval=len(no_bstar_nodes),
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                           progressbar.Percentage(), ' ', progressbar.Timer()])
    bar.start()
    default_num_reachable_nodes = {}
    for i, x in enumerate(no_bstar_nodes):
        T_A = bfs_hyperpath(x, set(), e, fstar)
        default_num_reachable_nodes[x] = len(T_A)
        bar.update(i)
    bar.finish()

    maxval = 10
    bar = progressbar.ProgressBar(maxval=maxval,
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                           progressbar.Percentage(), ' ', progressbar.Timer()])
    bar.start()
    str_to_write = ''
    for i in range(10):
        end = int(n * 0.03 * (i + 1))
        A = set(bc_sorted_idx[start:end])
        num_addi_reachable_nodes = 0
        for x in no_bstar_nodes:
            T_A = bfs_hyperpath(x, A, e, fstar)
            num_addi_reachable_nodes += (len(T_A) - default_num_reachable_nodes[x]) / len(no_bstar_nodes)
        str_to_write += str((i + 1) * 3) + '\t({0})\t'.format(len(A)) + str(num_addi_reachable_nodes) + '\n'
        bar.update(i)
    bar.finish()

    with fopen(reachpath, mode=mode) as f:
        if mode == 'w':
            f.write('k\t(#n)\t#addi. reachable nodes\n')
        f.write(str_to_write)

if __name__ == '__main__':
    usage = '''usage: python eval_reachability.py [filepath] [bc option] [epsilons] [etas]

    [bc option] = obc | bbc | lbbc | slbbc'''

    inCommand = sys.argv[1:]
    filepath = ''
    if len(inCommand) >= 1:
        filepath = inCommand[0]
    if len(inCommand) >= 2:
        if inCommand[1].lower() == 'obc':
            eval_reachability(filepath, bcoption='obc')
        elif inCommand[1].lower() == 'bbc':
            eval_reachability(filepath, bcoption='bbc')
        elif inCommand[1].lower() == 'lbbc':
            eval_reachability(filepath, bcoption='lbbc')
        elif inCommand[1].lower() == 'slbbc':
            if len(inCommand) >= 3:
                temp = inCommand[2].split(',')
                epsilons = []
                for eps in temp:
                    epsilons.append(float(eps))
                if len(inCommand) >= 4:
                    temp = inCommand[3].split(',')
                    etas = []
                    for eta in etas:
                        etas.append(float(eta))
                    eval_reachability(filepath, bcoption='slbbc', epsilons=epsilons, etas=etas)
                else:
                    eval_reachability(filepath, bcoption='slbbc', epsilons=epsilons)
            else:
                eval_reachability(filepath, bcoption='slbbc')
    else:
        print(usage)
