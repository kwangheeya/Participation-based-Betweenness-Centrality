from common import get_function_name,debug_print,read_hypergraph,read_bc,lcc_size


def coda_analysis():
    n, e = read_hypergraph("data/coda2.txt")
    bc = read_bc("bc/bch_apprx_coda2_0.5.txt", n)
    '''
    herbal_comp = read_seedset("data/coda2_herbal_compound_idmap.txt")
    for i in herbal_comp:
        if i < n:
            print(i, "\t", bc[i])
    '''
    idmap = read_idmap("data/coda2_idmap.txt")
    find_topAssociatedGene(n, bc, idmap,"data/coda2_asso_gene_idmap.txt")
    find_topGeneOrBody(n, bc, idmap,1,"data/coda2_asso_gene_idmap.txt")
    find_topAssociatedGene(n, bc, idmap,"data/coda2_diseaseGene_idmap.txt")
    find_topGeneOrBody(n, bc, idmap,1,"data/coda2_diseaseGene_idmap.txt")

def find_topAssociatedGene(n, bc, idmap, filename):
    asso_gene = read_seedset(filename)
    genebc = {}
    max_bc = 0
    for i in asso_gene:
        if i < n:
            genebc[i] = bc[i]
            if bc[i] > max_bc:
                max_bc = bc[i]
    print("max_bc", max_bc)
    i = 0
    for w in sorted(genebc, key=genebc.get, reverse=True):
        if i < 10:
            print(idmap[w], '\t',genebc[w])
        else:
            break
        i += 1

#"data/coda2_asso_gene_idmap.txt"
def find_topGeneOrBody(n, bc, idmap,option,filename):
    asso_gene = read_seedset(filename)
    genebc = {}
    max_bc = 0
    for i in asso_gene:
        if i < n:
            if bc[i] > 0:
                entity_name = idmap[i].split(' ')[option] #option = 0:gene / = 1:body
                if(entity_name in genebc):
                    genebc[entity_name].add(i)
                else:
                    genebc[entity_name] = set([i])
    averagebc= {}
    for key in genebc:
        sumbc = 0
        for eid in genebc[key]:
            sumbc += bc[eid]
        averagebc[key] = sumbc/ float(len(genebc[key]))


    i = 0
    for w in sorted(averagebc, key=averagebc.get, reverse=True):
        if i < 10:
            print(w, '\t',averagebc[w])
        else:
            break
        i += 1


def read_idmap(_filename):
    HYPERGRAPH_file = './' + str(_filename)
    try:
        HYPERGRAPH_fp = open(HYPERGRAPH_file, 'r');
    except IOError:
        print('Error: can\'t find ', HYPERGRAPH_file, ' or read data')
        return None
    else:
        print('Reading data', '[' + HYPERGRAPH_file + ']')

    idmap = {}
    line = HYPERGRAPH_fp.readline()

    for line in HYPERGRAPH_fp:
        temp = line[:-1].split('\t')
        idmap[int(temp[0])] = temp[1]

    print('> Done. #id map:', len(idmap))
    HYPERGRAPH_fp.close()
    return idmap

def read_seedset(_filename):
    HYPERGRAPH_file = './' + str(_filename)
    try:
        HYPERGRAPH_fp = open(HYPERGRAPH_file, 'r');
    except IOError:
        print('Error: can\'t find ', HYPERGRAPH_file, ' or read data')
        return None
    else:
        print('Reading data', '[' + HYPERGRAPH_file + ']')

    seed_set = set()

    for line in HYPERGRAPH_fp:
        temp = line[:-1].split('\t')
        seed_id = temp[0].replace("[", "")
        seed_id = seed_id.replace("]", "")
        seed_set.add(int(seed_id))

    print('> Done. #seed nodes:', len(seed_set))
    HYPERGRAPH_fp.close()
    return seed_set


if __name__ == '__main__':
	coda_analysis()