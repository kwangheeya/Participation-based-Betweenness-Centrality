import sys
from pbc import OBC, PBC


def get_bcpath(filepath, bcoption, eps=0.01, is_get_bc_hyperedge=False):
    filename = filepath.split('/')[-1]
    bc_path = 'result/bc/'
    if bcoption.lower() == 'slpbc':
        s_filename = filename.split('.')
        if is_get_bc_hyperedge:
            bc_path += 'h_SLPBC_' + s_filename[0] + '_' + str(eps).split('.')[-1] + '.' + s_filename[1]
        else:
            bc_path += 'SLPBC_' + s_filename[0] + '_' + str(eps).split('.')[-1] + '.' + s_filename[1]
    elif bcoption.lower() == 'naiveslpbc':
        s_filename = filename.split('.')
        if is_get_bc_hyperedge:
            bc_path += 'h_NAIVE_SLPBC_' + s_filename[0] + '_' + str(eps).split('.')[-1] + '.' + s_filename[1]
        else:
            bc_path += 'NAIVE_SLPBC_' + s_filename[0] + '_' + str(eps).split('.')[-1] + '.' + s_filename[1]
    elif bcoption.lower() == 'bbc':
        s_filename = filename.split('.')
        bc_path += 'BBC_' + s_filename[0] + '_' + str(eps).split('.')[-1]+'.'+s_filename[1]
    else:
        bc_path += bcoption.upper() + '_' + filename
    return bc_path


def compute_obc(filepath):
    obc = OBC()
    obc.load_graph(filepath)
    obc.compute()
    bc_path = get_bcpath(filepath, 'OBC')
    obc.store_bc(bc_path)


def compute_pbc(filepath, bcoption='pbc', epsilons=[0.01], etas=[0.1]):
    #pbc = PBC(bcoption=bcoption, is_get_bc_hyperedge=True)
    pbc = PBC(bcoption=bcoption)
    pbc.load_graph(filepath)
    if bcoption == 'slpbc':
        for eps in epsilons:
            for eta in etas:
                pbc.compute(eps=eps, eta=eta)
                bc_path = get_bcpath(filepath, bcoption, eps=eps)
                pbc.store_bc(bc_path)
                if pbc.is_get_bc_hyperedge:
                    bc_path = get_bcpath(filepath, bcoption, eps=eps, is_get_bc_hyperedge=pbc.is_get_bc_hyperedge)
                    pbc.store_bcE(bc_path)
    elif bcoption == 'naiveslpbc':
        for eps in epsilons:
            for eta in etas:
                pbc.compute(eps=eps, eta=eta)
                bc_path = get_bcpath(filepath, bcoption, eps=eps)
                pbc.store_bc(bc_path)
                if pbc.is_get_bc_hyperedge:
                    bc_path = get_bcpath(filepath, bcoption, eps=eps, is_get_bc_hyperedge=pbc.is_get_bc_hyperedge)
                    pbc.store_bcE(bc_path)
    elif bcoption == 'bbc':
        for eps in epsilons:
            pbc.compute(eps=eps)
            bc_path = get_bcpath(filepath, bcoption, eps=eps)
            pbc.store_bc(bc_path)
    else:
        pbc.compute()
        bc_path = get_bcpath(filepath, bcoption)
        pbc.store_bc(bc_path)


if __name__ == '__main__':
    usage = '''usage: python compute_bc.py [filepath] [bc option] [epsilons] [etas]
    
    [bc option] = obc | pbc | lpbc | slpbc | bbc | naiveslpbc
    
    Examples: 
        python compute_bc.py data/hlmn/hlmn.txt obc
        python compute_bc.py data/hlmn/hlmn.txt pbc
        python compute_bc.py data/hlmn/hlmn.txt lpbc
        python compute_bc.py data/hlmn/hlmn.txt slpbc 0.01
        python compute_bc.py data/hlmn/hlmn.txt slpbc 0.01,0.005
        python compute_bc.py data/hlmn/hlmn.txt slpbc 0.01 0.1'''

    inCommand = sys.argv[1:]
    filepath = ''
    if len(inCommand) >= 1:
        filepath = inCommand[0]

    if len(inCommand) >= 2:
        if inCommand[1].lower() == 'obc':
            compute_obc(filepath)
        elif inCommand[1].lower() == 'pbc':
            compute_pbc(filepath, bcoption='pbc')
        elif inCommand[1].lower() == 'lpbc':
            compute_pbc(filepath, bcoption='lpbc')
        elif inCommand[1].lower() == 'slpbc':
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
                    compute_pbc(filepath, bcoption='slpbc', epsilons=epsilons, etas=etas)
                else:
                    compute_pbc(filepath, bcoption='slpbc', epsilons=epsilons)
            else:
                compute_pbc(filepath, bcoption='slpbc')
        elif inCommand[1].lower() == 'naiveslpbc':
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
                    compute_pbc(filepath, bcoption='naiveslpbc', epsilons=epsilons, etas=etas)
                else:
                    compute_pbc(filepath, bcoption='naiveslpbc', epsilons=epsilons)
            else:
                compute_pbc(filepath, bcoption='naiveslpbc')
        elif inCommand[1].lower() == 'bbc':
            if len(inCommand) >= 3:
                temp = inCommand[2].split(',')
                epsilons = []
                for eps in temp:
                    epsilons.append(float(eps))
                compute_pbc(filepath, bcoption='bbc', epsilons=epsilons)
            else:
                compute_pbc(filepath, bcoption='bbc', epsilons=[0.5])
        elif inCommand[1].lower() == 'hpbc':
            compute_pbc(filepath, bcoption='hpbc')

    else:
        print(usage)
