import sys
from bbc import OBC, BBC, fopen
from compute_bc import get_bcpath
import progressbar


def get_attackpath(filepath, bcoption, eps=0.01):
    filename = filepath.split('/')[-1]
    _path = 'result/attack/v_'
    if bcoption.lower() == 'slbbc':
        s_filename = filename.split('.')
        _path += 'SLBBC_' + s_filename[0] + '_' + str(eps).split('.')[-1]+'.'+s_filename[1]
    else:
        _path += bcoption.upper() + '_' + filename
    return _path


def eval_robustness(filepath, bcoption='slbbc', epsilons=[0.01], etas=[0.1]):
    bc = None
    if bcoption == 'obc':
        bc = OBC()
    else:
        bc = BBC(bcoption=bcoption)

    bc.load_graph(filepath)

    if bcoption == 'slbbc':
        for eps in epsilons:
            for eta in etas:
                bc_path = get_bcpath(filepath, bcoption, eps=eps)
                bc.load_bc(bc_path)
                attackpath = get_attackpath(filepath, bcoption, eps=eps)
                print('+ Start attack robustness in {0} with eps={1} eta={2}'.format(bcoption.upper(), eps, eta))
                _attack_robustness(bc, attackpath, eps=eps, eta=eta)
    else:
        bc_path = get_bcpath(filepath, bcoption)
        bc.load_bc(bc_path)
        attackpath = get_attackpath(filepath, bcoption)
        print('+ Start attack robustness in {0}'.format(bcoption.upper()))
        _attack_robustness(bc, attackpath)


def _attack_robustness(bc, attackpath, eps=0.01, eta=0.1):
    n = bc.graph.n

    with fopen(attackpath, mode='a+') as f:
        removed_nodes = []
        f.seek(0)
        for line in f:
            nid_str = line.split('\t')[2]
            if nid_str.isdigit():
                removed_nodes.append(int(nid_str))

        possible_nodes = []
        if len(removed_nodes) > 0:
            for nid in removed_nodes[:-1]:
                bc.graph.remove_node(nid)
            bc.bc = None
            inodes, possible_nodes = delete_node(bc, removed_nodes[-1])
            print('\t>Initially removed nodes:', len(removed_nodes))
        if f.tell() == 0:
            f.write('lcc\tinodes\tindex\tbc_val\n')

        maxval = n-len(removed_nodes)
        bar = progressbar.ProgressBar(max_value=maxval, widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                            progressbar.Percentage(), ' ', progressbar.AdaptiveETA()])
        bar.start()
        for i in range(n-len(removed_nodes)):
            if not bc.bc:
                bc.compute(print_op=False, possible_nodes=possible_nodes, eps=eps, eta=eta)
            nid = bc.sorted_index()[0]
            bc_val = bc.bc[nid]
            inodes, possible_nodes = delete_node(bc, nid)
            lcc = bc.graph.lcc()
            str_to_write = '{0}\t{1}\t{2}\t{3}\n'.format(lcc, inodes, nid, bc_val)
            f.write(str_to_write)
            bar.update(i)
            if lcc == 1 or bc_val == 0:
                break
            else:
                bc.bc = None

        bar.finish()


def delete_node(bc, nid):
    bc.graph.remove_node(nid)
    e = bc.graph.e
    n = bc.graph.n
    non_inodes = bc.graph.find_nodes_having_edges()
    n_inodes = n - len(non_inodes)

    return n_inodes, list(non_inodes)


if __name__ == '__main__':
    usage = '''usage: python eval_robustness.py [filepath] [bc option] [epsilons] [etas]

    [bc option] = obc | bbc | lbbc | slbbc'''

    inCommand = sys.argv[1:]
    filepath = ''
    if len(inCommand) >= 1:
        filepath = inCommand[0]
    if len(inCommand) >= 2:
        if inCommand[1].lower() == 'obc':
            eval_robustness(filepath, bcoption='obc')
        elif inCommand[1].lower() == 'bbc':
            eval_robustness(filepath, bcoption='bbc')
        elif inCommand[1].lower() == 'lbbc':
            eval_robustness(filepath, bcoption='lbbc')
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
                    eval_robustness(filepath, bcoption='slbbc', epsilons=epsilons, etas=etas)
                else:
                    eval_robustness(filepath, bcoption='slbbc', epsilons=epsilons)
            else:
                eval_robustness(filepath, bcoption='slbbc')
    else:
        print(usage)
