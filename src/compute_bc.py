import sys
from bbc import Bhypergraph, Graph

def compute_obc(filepath):
    graph = Graph(filepath)
    print(graph.lcc())

def compute_bbc(filepath, bcoption = 'e', epsilons=[0.01], etas=[0.1]):
    bhyper = Bhypergraph(filepath)
    print(bhyper.lcc())

if __name__ == '__main__':
    usage = 'python compute_bc.py [filepath] [bc option] [epsilons] [etas]'

    inCommand = sys.argv[1:]
    filepath = ''
    if len(inCommand) >= 1:
        filepath = inCommand[0]

    if len(inCommand) >= 2:
        if inCommand[1] == 'obc':
            compute_obc(filepath)
        elif inCommand[1] == 'bbc':
            compute_bbc(filepath, bcoption='e')
        elif inCommand[1] == 'lbbc':
            compute_bbc(filepath, bcoption='l')
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
                    compute_bbc(filepath, bcoption='s', epsilons=epsilons, etas=etas)
                else:
                    compute_bbc(filepath, bcoption='s', epsilons=epsilons)
            else:
                compute_bbc(filepath, bcoption='s')
    else:
        print(usage)
