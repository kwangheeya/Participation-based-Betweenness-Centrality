import sys
from bbc import Bhypergraph, Graph, OBC


def get_bcpath(filepath, bcoption, compute_opt='s', eps=0.01):
    filename = filepath.split('/')[-1]
    bc_path = 'result/bc/'
    if bcoption.lower() == 'slbbc':
        s_filename = filename.split('.')
        bc_path += 'slBBC_' + s_filename[0] + '_' + str(eps).split('.')[-1]+'.'+s_filename[1]
    else:
        bc_path += bcoption + '_' + filename
    return bc_path


def compute_obc(filepath):
    obc = OBC()
    obc.load_graph(filepath)
    obc.compute()
    bc_path = get_bcpath(filepath, 'OBC')
    obc.store_bc(bc_path)


def compute_bbc(filepath, compute_opt='e', epsilons=[0.01], etas=[0.1]):
    bhyper = Bhypergraph(filepath)


if __name__ == '__main__':
    usage = '''usage: python compute_bc.py [filepath] [bc option] [epsilons] [etas]
    
    [bc option] = obc | bbc | lbbc | slbbc
    
    Example: python compute_bc.py data/hlmn/hlmn.txt obc
    Example: python compute_bc.py data/hlmn/hlmn.txt slbbc 0.01
    Example: python compute_bc.py data/hlmn/hlmn.txt slbbc 0.01,0.005
    Example: python compute_bc.py data/hlmn/hlmn.txt slbbc 0.01 0.1'''

    inCommand = sys.argv[1:]
    filepath = ''
    if len(inCommand) >= 1:
        filepath = inCommand[0]

    if len(inCommand) >= 2:
        if inCommand[1] == 'obc':
            compute_obc(filepath)
        elif inCommand[1] == 'bbc':
            compute_bbc(filepath, compute_opt='e')
        elif inCommand[1] == 'lbbc':
            compute_bbc(filepath, compute_opt='l')
        elif inCommand[1] == 'slbbc':
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
                    compute_bbc(filepath, compute_opt='s', epsilons=epsilons, etas=etas)
                else:
                    compute_bbc(filepath, compute_opt='s', epsilons=epsilons)
            else:
                compute_bbc(filepath, compute_opt='s')
    else:
        print(usage)
