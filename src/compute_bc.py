import sys
from bbc import OBC, BBC


def get_bcpath(filepath, bcoption, eps=0.01):
    filename = filepath.split('/')[-1]
    bc_path = 'result/bc/'
    if bcoption.lower() == 'slbbc':
        s_filename = filename.split('.')
        bc_path += 'SLBBC_' + s_filename[0] + '_' + str(eps).split('.')[-1]+'.'+s_filename[1]
    else:
        bc_path += bcoption.upper() + '_' + filename
    return bc_path


def compute_obc(filepath):
    obc = OBC()
    obc.load_graph(filepath)
    obc.compute()
    bc_path = get_bcpath(filepath, 'OBC')
    obc.store_bc(bc_path)


def compute_bbc(filepath, bcoption='bbc', epsilons=[0.01], etas=[0.1]):
    bbc = BBC(bcoption=bcoption)
    bbc.load_graph(filepath)
    if bcoption == 'slbbc':
        for eps in epsilons:
            for eta in etas:
                bbc.compute(eps=eps, eta=eta)
                bc_path = get_bcpath(filepath, bcoption, eps=eps)
                bbc.store_bc(bc_path)
    else:
        bbc.compute()
        bc_path = get_bcpath(filepath, bcoption)
        bbc.store_bc(bc_path)


if __name__ == '__main__':
    usage = '''usage: python compute_bc.py [filepath] [bc option] [epsilons] [etas]
    
    [bc option] = obc | bbc | lbbc | slbbc
    
    Examples: 
        python compute_bc.py data/hlmn/hlmn.txt obc
        python compute_bc.py data/hlmn/hlmn.txt slbbc 0.01
        python compute_bc.py data/hlmn/hlmn.txt slbbc 0.01,0.005
        python compute_bc.py data/hlmn/hlmn.txt slbbc 0.01 0.1'''

    inCommand = sys.argv[1:]
    filepath = ''
    if len(inCommand) >= 1:
        filepath = inCommand[0]

    if len(inCommand) >= 2:
        if inCommand[1].lower() == 'obc':
            compute_obc(filepath)
        elif inCommand[1].lower() == 'bbc':
            compute_bbc(filepath, bcoption='bbc')
        elif inCommand[1].lower() == 'lbbc':
            compute_bbc(filepath, bcoption='lbbc')
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
                    compute_bbc(filepath, bcoption='slbbc', epsilons=epsilons, etas=etas)
                else:
                    compute_bbc(filepath, bcoption='slbbc', epsilons=epsilons)
            else:
                compute_bbc(filepath, bcoption='slbbc')
    else:
        print(usage)
